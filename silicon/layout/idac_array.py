#!/usr/bin/env python3
"""FDVPD current DAC, unit-source array -- sky130A layout generator.

Stage 1 of the DAC layout: every unit current source of the DAC and of its
reference diode, as one matched array, with each unit's drain routed out to
its segment's net.  Switches and load resistors are the next stage.

What gets drawn is fixed by layout/idac_netlist.py: 1023 output units in 21
segments (15 x 64 thermometer, 32..1 binary), 64 fixed reference-diode units
and 127 trim units in 7 binary groups -- 1214 devices of W/L = 0.42/24.

The array is built from one idea: all units share a gate (ibias), so a single
24 um wide poly sheet per column with 0.42 um diffusion stripes running under
it IS a column of W/L = 0.42/24 transistors, one per stripe.  Every unit then
sees the same poly, the same neighbours and the same contact geometry.

    x ->   tap | S [ poly sheet, L = 24 um ] D | channel | D [ sheet ] S | tap
                          column 2k                         column 2k+1

  * Sources face outward and merge with the p-tap strips into one m1 vss
    plane, so the sources need no routing at all, and every device is within
    about 13 um of a tap.
  * Drains face a shared routing channel: one m2 track per net, a left half
    for the left column and a right half for the right one, so drain lines at
    the same height from the two columns can never meet.
  * Below the array every track drops onto an m3 bus per net; the buses are
    the ports.

Matching.  Units are handed out in point-symmetric PAIRS about the array
centre, so every segment's linear gradient cancels exactly.  The pairs are
dealt over a spread ordering of the array by largest remainder, so each
segment also samples the whole array rather than a patch of it.  The two
1-unit nets (binary bit 0 and trim bit 0) share the centre pair.  Two dummy
stripes close each column top and bottom.

    python3 idac_array.py              # writes idac_array.gds
    python3 idac_array.py --png x.png  # also renders it
"""
import argparse
import os
import sys

import klayout.db as kdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import idac_netlist as dac  # noqa: E402
from ramp_sink import (LAYERS, Drawer, _array, snap, render,  # noqa: E402
                       LIC, LIC_P, MCON_P, VIA, VIA2, SDM_ENC, RING_TAP,
                       GB_H)

CELL = "idac_array"

# --- geometry (nm) -----------------------------------------------------------
UW, UL = int(round(dac.UNIT_W * 1000)), int(round(dac.UNIT_L * 1000))  # unit W, L
PITCH = UW + 270         # stripe pitch: diffusion + diff spacing 270
SDX = 330               # diff past poly: licon.5c 60 + licon 170 + 100 to gate
LIC_X = 60 + LIC // 2   # licon centre from the diffusion end
M1H = 260               # drain m1 line: via1 150 + 2 x 55 (via.4a)
TRK_W = 320             # m2 track: via1 150 + 2 x 85
TRK_P = 500             # m2 track pitch
BUS_W = 330             # m3 bus: via2 200 + 2 x 65
BUS_P = 700
GATE_GAP = 250          # outermost diff -> gate-bar poly
TAP_GAP = 600           # diffusion -> tap strip
CHAN_MARGIN = 600       # drain-side tap -> first track
ENC_M1_MCON = 60        # met1.5: one direction
ENC_M1_VIA = 85         # via.5a: one direction

N_COLS = dac.ARRAY_COLS
EDGE_DUMMIES = dac.ARRAY_EDGE_DUMMIES


def assign(rows):
    """Map (col, row) -> net for the inner (non-edge) positions."""
    cols = N_COLS
    cx, cy = (cols - 1) / 2, (rows - 1) / 2

    def mirror(p):
        return (cols - 1 - p[0], rows - 1 - p[1])

    # one representative of every point-symmetric pair: the lower half by row,
    # and on the centre row (odd rows only) the left half by column
    half = [(c, r) for r in range(rows) for c in range(cols)
            if (r, c) < mirror((c, r))[::-1]]
    # spread order: bit-reversed row index, so consecutive picks land far apart
    bits = max(1, (rows - 1).bit_length())

    def rev(r):
        return int(f"{r:0{bits}b}"[::-1], 2)
    half.sort(key=lambda p: (rev(p[1]), (p[0] * 7) % cols))

    items = dac.array_units()
    odd = [n for n, k in items if k % 2]
    if len(odd) != 2:
        raise ValueError(f"expected exactly two odd nets, got {odd}")
    pairs = {n: k // 2 for n, k in items}
    total_pairs = sum(pairs.values()) + 1        # + the shared odd pair
    if total_pairs > len(half):
        raise ValueError(f"{total_pairs} pairs do not fit in {len(half)}")
    pairs["dummy"] = len(half) - total_pairs

    # the shared odd pair takes the representative nearest the centre
    centre = min(half, key=lambda p: (p[0] - cx) ** 2 + (p[1] - cy) ** 2)
    half.remove(centre)
    where = {centre: odd[0], mirror(centre): odd[1]}

    # largest-remainder dealing: at every step give the pair to the net that
    # is furthest behind its share -- each net spreads evenly over `half`
    n_total = len(half)
    given = {n: 0 for n in pairs}
    for i, p in enumerate(half):
        net = max(pairs, key=lambda n: pairs[n] * (i + 1) / n_total - given[n])
        given[net] += 1
        where[p] = net
        where[mirror(p)] = net
    assert all(given[n] == pairs[n] for n in pairs), (given, pairs)
    return where


def build(ly=None):
    """Draw the array.  Into ``ly`` if given (as a subcell of a bigger layout),
    else into a fresh layout.  Returns (layout, cell, stats, geometry), where
    geometry is what a parent needs to connect to it: the m3 bus heights, the
    column spans (the only x ranges with no m2 in the bus region), and the
    extent."""
    if ly is None:
        ly = kdb.Layout()
        ly.dbu = 0.001
    top = ly.create_cell(CELL)
    d = Drawer(ly, top)

    items = dac.array_units()
    n_units = sum(k for _, k in items)
    inner_rows, rows = dac.array_rows()
    where = assign(inner_rows)

    chan_nets = [n for n, _ in items] + ["vss"]
    half_w = len(chan_nets) * TRK_P
    dev_w = SDX + UL + SDX                        # one stripe, end to end
    col_h = (rows - 1) * PITCH + UW

    # ---- x plan -------------------------------------------------------------
    x = RING_TAP + TAP_GAP                        # left tap strip at x = 0
    cols = []                                      # (x_diff0, drains_right)
    chans = []                                     # (x_left_half0, x_right_half0)
    taps = [0]
    drain_taps = []                                # tap strips beside the drain ends
    for p in range(N_COLS // 2):
        cols.append((x, True))
        drain_taps.append(x + dev_w + TAP_GAP)
        x_ch = x + dev_w + TAP_GAP + RING_TAP + CHAN_MARGIN
        chans.append((x_ch, x_ch + half_w))
        x = x_ch + 2 * half_w + CHAN_MARGIN
        drain_taps.append(x)
        x += RING_TAP + TAP_GAP
        cols.append((x, False))
        x += dev_w + TAP_GAP
        taps.append(x)
        x += RING_TAP + TAP_GAP
    x_end = taps[-1] + RING_TAP

    y0 = 0                                         # bottom of stripe 0
    y_gate_bot = y0 - GATE_GAP - GB_H
    y_gate_top = y0 + col_h + GATE_GAP

    def track_x(col, net):
        ch = chans[col // 2]
        base = ch[0] if col % 2 == 0 else ch[1]
        return base + TRK_P / 2 + chan_nets.index(net) * TRK_P

    # ---- columns -----------------------------------------------------------
    for ci, (xd0, drains_right) in enumerate(cols):
        xd1 = xd0 + dev_w
        xp0, xp1 = xd0 + SDX, xd1 - SDX
        xs = xd0 + LIC_X if drains_right else xd1 - LIC_X     # source licon x
        xdr = xd1 - LIC_X if drains_right else xd0 + LIC_X    # drain licon x
        # poly sheet with gate bars top and bottom
        d.rect("poly", xp0, y_gate_bot, xp1, y_gate_top + GB_H)
        d.rect("nsdm", xd0 - SDM_ENC, y0 - SDM_ENC, xd1 + SDM_ENC,
               y0 + col_h + SDM_ENC)
        for yb in (y_gate_bot, y_gate_top):
            gy = yb + GB_H / 2
            d.rect("npc", xp0, yb, xp1, yb + GB_H)
            for gx in _array(xp0 + 100, xp1 - 100, LIC, LIC_P):
                d.square("licon", gx, gy, LIC)
            d.rect("li", xp0 + 20, yb + 20, xp1 - 20, yb + GB_H - 20)
            for gx in _array(xp0 + 300, xp1 - 300, LIC, MCON_P):
                d.square("mcon", gx, gy, LIC)
            # m1 gate strip, out to the ibias track in this column's channel
            tx = track_x(ci, "ibias")
            if drains_right:
                d.rect("m1", xp0 + 200, gy - M1H / 2, tx + VIA / 2 + ENC_M1_VIA, gy + M1H / 2)
            else:
                d.rect("m1", tx - VIA / 2 - ENC_M1_VIA, gy - M1H / 2, xp1 - 200, gy + M1H / 2)
            d.square("via", tx, gy, VIA)
        # source li bar + mcons; the m1 vss plane is drawn with the taps
        d.rect("li", xs - LIC / 2, y0 + UW / 2 - LIC / 2 - 80, xs + LIC / 2,
               y0 + col_h - UW / 2 + LIC / 2 + 80)
        for r in range(rows):
            yc = y0 + r * PITCH + UW / 2
            d.rect("diff", xd0, yc - UW / 2, xd1, yc + UW / 2)
            d.square("licon", xs, yc, LIC)
            d.square("mcon", xs, yc, LIC)
            # drain: licon, li pad, mcon, m1 line to the net's track
            inner = EDGE_DUMMIES <= r < rows - EDGE_DUMMIES
            net = where.get((ci, r - EDGE_DUMMIES), "dummy") if inner else "dummy"
            net = "vss" if net == "dummy" else net
            tx = track_x(ci, net)
            d.square("licon", xdr, yc, LIC)
            d.rect("li", xdr - LIC / 2 - 80, yc - LIC / 2,
                   xdr + LIC / 2 + 80, yc + LIC / 2)
            d.square("mcon", xdr, yc, LIC)
            if drains_right:
                d.rect("m1", xdr - LIC / 2 - ENC_M1_MCON, yc - M1H / 2,
                       tx + VIA / 2 + ENC_M1_VIA, yc + M1H / 2)
            else:
                d.rect("m1", tx - VIA / 2 - ENC_M1_VIA, yc - M1H / 2,
                       xdr + LIC / 2 + ENC_M1_MCON, yc + M1H / 2)
            d.square("via", tx, yc, VIA)

    # ---- taps + m1 vss planes (tap strip merged with the source columns) ----
    y_lo, y_hi = y_gate_bot - 400, y_gate_top + GB_H + 400
    bus_y = {n: y_lo - 1200 - i * BUS_P for i, n in enumerate(chan_nets)}
    for i, tx0 in enumerate(taps):
        tx1 = tx0 + RING_TAP
        d.rect("tap", tx0, y_lo, tx1, y_hi)
        d.rect("psdm", tx0 - SDM_ENC, y_lo - SDM_ENC, tx1 + SDM_ENC, y_hi + SDM_ENC)
        d.rect("li", tx0 + 40, y_lo + 40, tx1 - 40, y_hi - 40)
        for y in _array(y_lo + 120, y_hi - 120, LIC, LIC_P):
            d.square("licon", tx0 + RING_TAP / 2, y, LIC)
        for y in _array(y_lo + 120, y_hi - 120, LIC, MCON_P):
            d.square("mcon", tx0 + RING_TAP / 2, y, LIC)
        # plane: from the source contacts on either side, across the tap
        left = [c for c in cols if not c[1] and c[0] + dev_w < tx0]
        right = [c for c in cols if c[1] and c[0] > tx1]
        px0 = (left[-1][0] + dev_w - LIC_X - LIC / 2 - 30) if left else tx0
        px1 = (right[0][0] + LIC_X + LIC / 2 + 30) if right else tx1
        d.rect("m1", px0, y0 - 100, px1, y0 + col_h + 100)
        d.rect("m1", tx0, y_lo, tx1, y_hi)
        # down to the vss bus
        d.rect("m1", tx0, bus_y["vss"] - 200, tx1, y_lo)
        cx = tx0 + RING_TAP / 2
        d.square("via", cx, bus_y["vss"], VIA)
        d.square("m2", cx, bus_y["vss"], 400)
        d.square("via2", cx, bus_y["vss"], VIA2)

    # ---- drain-side taps.  LU.2 wants every n-diffusion within 15 um of a
    # p-tap, and the drain ends sit 25 um from the outer ones.  The drain m1
    # lines cross these strips, so they are strapped in li alone and reach vss
    # below the array, where m1 is free.
    for tx0 in drain_taps:
        tx1 = tx0 + RING_TAP
        cx = tx0 + RING_TAP / 2
        d.rect("tap", tx0, y_lo, tx1, y_hi)
        d.rect("psdm", tx0 - SDM_ENC, y_lo - SDM_ENC, tx1 + SDM_ENC, y_hi + SDM_ENC)
        for y in _array(y_lo + 120, y_hi - 120, LIC, LIC_P):
            d.square("licon", cx, y, LIC)
        d.rect("li", tx0 + 40, bus_y["vss"] - 200, tx1 - 40, y_hi - 40)
        d.square("mcon", cx, bus_y["vss"], LIC)
        d.rect("m1", cx - 200, bus_y["vss"] - 200, cx + 200, bus_y["vss"] + 200)
        d.square("via", cx, bus_y["vss"], VIA)
        d.square("m2", cx, bus_y["vss"], 400)
        d.square("via2", cx, bus_y["vss"], VIA2)

    # ---- channel tracks (m2) down to their buses (m3) -----------------------
    y_top_trk = y_gate_top + GB_H / 2 + 200
    for ci in range(N_COLS):
        for net in chan_nets:
            tx = track_x(ci, net)
            by = bus_y[net]
            d.rect("m2", tx - TRK_W / 2, by - 200, tx + TRK_W / 2, y_top_trk)
            d.square("via2", tx, by, VIA2)

    x_bus0 = -2000
    for net in chan_nets:
        by = bus_y[net]
        d.rect("m3", x_bus0, by - BUS_W / 2, x_end, by + BUS_W / 2)
        d.rect("m3pin", x_bus0, by - BUS_W / 2, x_bus0 + 600, by + BUS_W / 2)
        d.label("m3lbl", net, x_bus0 + 300, by)

    full = d.bbox(list(LAYERS))
    d.rect("bound", full.left, full.bottom, full.right, full.top)
    stats = {"rows": rows, "inner_rows": inner_rows, "units": n_units,
             "dummies": N_COLS * rows - n_units}
    geo = {"bus_y": bus_y, "bus_w": BUS_W, "x0": x_bus0, "x1": x_end,
           "cols": [(c[0], c[0] + dev_w) for c in cols],
           "bottom": full.bottom, "top": full.top}
    return ly, top, stats, geo


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("-o", "--out", default=os.path.join(here, CELL + ".gds"))
    ap.add_argument("--png", help="also render the layout to this PNG")
    args = ap.parse_args()
    ly, top, stats, _ = build()
    ly.write(args.out)
    bb = top.dbbox()
    print(f"wrote {args.out}: {bb.width():.1f} x {bb.height():.1f} um, "
          f"{N_COLS} x {stats['rows']} stripes, {stats['units']} units, "
          f"{stats['dummies']} dummies")
    if args.png:
        render(args.out, args.png)


if __name__ == "__main__":
    main()
