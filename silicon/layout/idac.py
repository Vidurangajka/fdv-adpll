#!/usr/bin/env python3
"""FDVPD current DAC -- top cell, sky130A layout generator.

Stage 2 of the DAC: the verified unit array (idac_array.py) as a subcell, and
below it a strip holding everything else the netlist (idac_netlist.py) has:

    21 cascodes       one per segment, right above its switch pair
    21 switch pairs   steer each segment to outp or outn
     7 trim switches  select the reference-diode trim units
     2 R_D            res_high_po_0p35, drawn by magic's own PDK generator

                 idac_array  (m3 buses along its bottom edge)
      | | | |    m2 drops, one per bus, straight down under the unit columns
      ======     outn bus (m4)
      ======     outp bus (m3) | ibias bus (m3), side by side
      ------  p-tap ring -----------------------------------------------
               [ cascode x21 ........................ ]  bottoms aligned
               ======  cbias (m1) through every gate bar
               ======  p-tap strip (LU.2: cascodes are up to 13 um tall)
    [trim x7]  [ pair x21 ........................... ]   [R_D x2]
      | | | |    gates in m2, down to pins on the bottom edge

Each cascode and switch pair sits directly under its segment's drop, so the
source nets need no horizontal routing at all.  The cascode has to be HERE,
next to the switch: the routed source net is 170-250 fF, and without a
cascode it is the node that moves when the switch changes sides (post-layout
DNL 0.42 LSB; 0.01 with).  A pair's two drains climb past the cascode on
different metals -- m3 to outp, m4 to outn -- which is why outn is on m4.
The drops run only under the unit columns: that is the one x range where the
array's bus region has no m2 of its own (idac_array.build returns the spans).

    python3 idac.py -o build/idac/idac.gds --png build/idac/idac.png
"""
import argparse
import os
import subprocess
import sys

import klayout.db as kdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import idac_array  # noqa: E402
import idac_netlist as dac  # noqa: E402
from ramp_sink import (LAYERS, Drawer, render, ring,  # noqa: E402
                       LIC, LIC_P, MCON_P, VIA, VIA2, SDM_ENC, RING_TAP,
                       _array)

LAYERS.update({"m2lbl": (69, 5), "m2pin": (69, 16), "via3": (70, 44),
               "m4": (71, 20), "m4lbl": (71, 5), "m4pin": (71, 16)})

CELL = "idac"

# --- device geometry (nm), shared by the switch pair and the trim switch ----
W, L = 420, 150            # every switch is 0.42 / 0.15
SD = 360                   # S/D column: licon 170 + 95 to each gate
HEAD_W, HEAD_H = 290, 370  # poly gate head: licon 170 + 60 / + 100
HEAD_GAP = 250             # diffusion -> head
M1P_W, M1P_H = 260, 400    # m1 pad: via1 150 + 2 x 55 across, 2 x 125 along
M2P_W, M2P_H = 320, 370    # m2 pad: via1 / via2 enclosure
M3_W = 330                 # m3 / m4 line: via 200 + 2 x 65
M2_W = 320
VIA3 = 200

PITCH = 3000               # drop / cell pitch along the strip
BUS_H = 400                # outp / outn / ibias / vdd bus height


def col_x(i):
    """Centre of S/D column i in a row that starts at x = 0."""
    return i * (SD + L) + SD / 2


def device_row(d, x0, y0, n_gates, w=W):
    """One diffusion of width ``w`` with ``n_gates`` fingers; gate heads below.

    Returns (column centres, gate centres, head centre y, via y).  The S/D
    columns get their contacts and a via1 at mid-height; what the via goes up
    to is left to the caller.
    """
    n_cols = n_gates + 1
    x1 = x0 + n_cols * SD + n_gates * L
    d.rect("diff", x0, y0, x1, y0 + w)
    d.rect("nsdm", x0 - SDM_ENC, y0 - SDM_ENC, x1 + SDM_ENC, y0 + w + SDM_ENC)
    cols = [x0 + col_x(i) for i in range(n_cols)]
    gates = [x0 + (i + 1) * SD + i * L + L / 2 for i in range(n_gates)]
    hb = y0 - HEAD_GAP - HEAD_H
    hy = hb + HEAD_H / 2
    for gx in gates:
        d.rect("poly", gx - L / 2, hb, gx + L / 2, y0 + w + 130)
        d.rect("poly", gx - HEAD_W / 2, hb, gx + HEAD_W / 2, hb + HEAD_H)
        d.square("licon", gx, hy, LIC)
        d.rect("li", gx - LIC / 2, hy - LIC / 2 - 80, gx + LIC / 2, hy + LIC / 2 + 80)
        d.square("mcon", gx, hy, LIC)
        d.rect("m1", gx - M1P_W / 2, hy - M1P_H / 2, gx + M1P_W / 2, hy + M1P_H / 2)
        d.square("via", gx, hy, VIA)
    d.rect("npc", gates[0] - HEAD_W / 2 - 40, hb, gates[-1] + HEAD_W / 2 + 40,
           hb + HEAD_H)
    yc = y0 + w / 2
    ys = _array(y0 + 60, y0 + w - 60, LIC, LIC_P)
    li0, li1 = ys[0] - LIC / 2 - 80, ys[-1] + LIC / 2 + 80
    for cx in cols:
        for y in ys:
            d.square("licon", cx, y, LIC)
        d.rect("li", cx - LIC / 2, li0, cx + LIC / 2, li1)
        for y in _array(li0, li1, LIC, MCON_P):
            d.square("mcon", cx, y, LIC)
        d.rect("m1", cx - M1P_W / 2, min(li0 - 60, yc - M1P_H / 2),
               cx + M1P_W / 2, max(li1 + 60, yc + M1P_H / 2))
        d.square("via", cx, yc, VIA)
    return cols, gates, hy, yc


def up_m2(d, x, y, y_top):
    """m2 line from a column's via1 straight up to y_top."""
    d.rect("m2", x - M2_W / 2, y - M2P_H / 2, x + M2_W / 2, y_top)


def up_m3(d, x, y, y_top):
    """via1 -> m2 pad -> via2 -> m3 line up to y_top."""
    d.rect("m2", x - M2P_W / 2, y - M2P_H / 2, x + M2P_W / 2, y + M2P_H / 2)
    d.square("via2", x, y, VIA2)
    d.rect("m3", x - M3_W / 2, y - 200, x + M3_W / 2, y_top)


def up_m4(d, x, y, y_top):
    """via1 -> m2 -> via2 -> m3 pad -> via3 -> m4 line up to y_top."""
    d.rect("m2", x - M2P_W / 2, y - M2P_H / 2, x + M2P_W / 2, y + M2P_H / 2)
    d.square("via2", x, y, VIA2)
    d.rect("m3", x - M3_W / 2, y - 400, x + M3_W / 2, y + 400)   # met3.6 area
    d.square("via3", x, y, VIA3)
    d.rect("m4", x - M3_W / 2, y - 200, x + M3_W / 2, y_top)


def gate_down(d, x, hy, y_pin, net):
    """Gate head via1 -> m2 line down to a pin on the bottom edge."""
    d.rect("m2", x - M2_W / 2, y_pin, x + M2_W / 2, hy + M1P_H / 2)
    d.rect("m2pin", x - M2_W / 2, y_pin, x + M2_W / 2, y_pin + 500)
    d.label("m2lbl", net, x, y_pin + 250)


CL = 500                   # cascode L


def cascode(d, xs, yb, wf):
    """Two-finger cascode [c | cbias | x | cbias | c], centred on x = xs.

    The outer columns (source, the routed array net) are strapped together in
    m1 above the diffusion, and the strap takes via1 at xs for the drop to come
    down onto.  The middle column (drain) takes via1 at its bottom, for an m2
    straight down to the switch pair below.  The gate bar sits under the
    diffusion, centred on the row's cbias line, which the caller draws.
    Returns (drain via y, strap y).
    """
    x0 = xs - (SD + CL + SD / 2)
    x1 = x0 + 3 * SD + 2 * CL
    d.rect("diff", x0, yb, x1, yb + wf)
    d.rect("nsdm", x0 - SDM_ENC, yb - SDM_ENC, x1 + SDM_ENC, yb + wf + SDM_ENC)
    cols = [x0 + i * (SD + CL) + SD / 2 for i in range(3)]
    gates = [x0 + (i + 1) * SD + i * CL + CL / 2 for i in range(2)]
    bar_top = yb - HEAD_GAP
    bar_bot = bar_top - HEAD_H
    for gx in gates:
        d.rect("poly", gx - CL / 2, bar_bot, gx + CL / 2, yb + wf + 130)
    bx0, bx1 = gates[0] - CL / 2, gates[1] + CL / 2
    d.rect("poly", bx0, bar_bot, bx1, bar_top)
    d.rect("npc", bx0, bar_bot, bx1, bar_top)
    by = (bar_bot + bar_top) / 2
    for x in _array(bx0 + 100, bx1 - 100, LIC, LIC_P):
        d.square("licon", x, by, LIC)
    d.rect("li", bx0 + 20, bar_bot + 20, bx1 - 20, bar_top - 20)
    for x in _array(bx0 + 150, bx1 - 150, LIC, MCON_P):
        d.square("mcon", x, by, LIC)
    ys = _array(yb + 60, yb + wf - 60, LIC, LIC_P)
    li0, li1 = ys[0] - LIC / 2 - 80, ys[-1] + LIC / 2 + 80
    y_strap = yb + wf + 400
    for i, cx in enumerate(cols):
        for y in ys:
            d.square("licon", cx, y, LIC)
        d.rect("li", cx - LIC / 2, li0, cx + LIC / 2, li1)
        for y in _array(li0, li1, LIC, MCON_P):
            d.square("mcon", cx, y, LIC)
        # met1.5: 60 past the end mcons along the column
        top = y_strap + M1H_S / 2 if i != 1 else li1 + 60
        d.rect("m1", cx - M1P_W / 2, li0 - 60, cx + M1P_W / 2, top)
    d.rect("m1", cols[0] - M1P_W / 2, y_strap - M1H_S / 2, cols[2] + M1P_W / 2,
           y_strap + M1H_S / 2)
    d.square("via", xs, y_strap, VIA)
    d.square("via", xs, ys[0], VIA)
    return ys[0], y_strap


M1H_S = 260                # m1 strap height


def make_rd(out_dir):
    """Both R_D from magic's own generator, as one two-element array.

    A lone 0.35 um high-sheet resistor is not DRC-legal as generated: its rpm
    marker comes out 0.75 um wide against rpm.1's 1.27 um, even with its guard
    ring.  Two side by side share one marker wide enough, one guard ring, and
    are the matched pair R_P / R_N should be anyway.
    """
    gds = os.path.join(out_dir, "idac_rd.gds")
    tcl = os.path.join(out_dir, "idac_rd.tcl")
    with open(tcl, "w") as fh:
        fh.write(f"""load idac_rd
box 0 0 0 0
set p [dict merge [sky130::{dac.RES}_defaults] \
{{w {dac.RD_W} l {dac.RD_L} nx 2}}]
sky130::{dac.RES}_draw $p
gds write {os.path.basename(gds)}
quit -noprompt
""")
    rc = os.path.join(os.environ["PDK_ROOT"], "sky130A", "libs.tech", "magic",
                      "sky130A.magicrc")
    subprocess.run(["magic", "-dnull", "-noconsole", "-rcfile", rc,
                    os.path.basename(tcl)], cwd=out_dir, check=True,
                   capture_output=True)
    return gds


def build(out_dir):
    ly = kdb.Layout()
    ly.dbu = 0.001
    top = ly.create_cell(CELL)
    d = Drawer(ly, top)

    # ---- the array, and where its buses are ------------------------------
    _, arr, _, geo = idac_array.build(ly)
    top.insert(kdb.CellInstArray(arr.cell_index(), kdb.Trans()))
    bus_y = geo["bus_y"]

    # ---- the resistor cell ----------------------------------------------
    rd_ly = kdb.Layout()
    rd_ly.read(make_rd(out_dir))
    rd_src = rd_ly.cell("idac_rd")
    rd = ly.create_cell("idac_rd")
    rd.copy_tree(rd_src)
    # geometry read from the generated cell, whatever resistor it is: the
    # terminal labels sit on the terminal strips, B on the guard ring
    lbl = {}
    for li in rd_ly.layer_indexes():
        for sh in rd_src.shapes(li).each():
            if sh.is_text():
                lbl[sh.text_string] = (sh.text_pos.x, sh.text_pos.y)
    bb = rd_src.bbox()
    RD_X = abs(lbl["R1_1"][0])         # the two resistors' x, +- from centre
    RD_B = abs(lbl["B"][1])            # guard ring bottom side, below centre
    RD_HW, RD_HH = bb.right, bb.top    # half width / height incl. guard ring
    # terminal pad y: the middle of the m1 strip under the R1_1 label
    m1 = rd_ly.find_layer(68, 20)
    pt = kdb.Point(*lbl["R1_1"])
    strips = [sh.bbox() for sh in rd_src.shapes(m1).each() if sh.bbox().contains(pt)]
    RD_TERM = strips[0].center().y if strips else abs(lbl["R1_1"][1])

    # ---- drop positions: only under the unit columns ----------------------
    slots = []
    for c0, c1 in geo["cols"]:
        x = c0 + 1500
        while x <= c1 - 1500:
            slots.append(x)
            x += PITCH
    trims = [f"dt{k}" for k in range(len(dac.REF_TRIM))]
    segs = [n for n, _ in dac.segments()]
    oss = [n for n, _ in dac.os_segments()]
    need = len(trims) + 1 + len(segs) + len(oss) + 3  # + ibias drop + vss drops
    if len(slots) < need:
        raise ValueError(f"{need} drops, only {len(slots)} slots under the columns")
    slot = iter(slots)
    vss_x = [next(slot) for _ in range(3)]         # first: inside the ring for sure
    trim_x = [next(slot) for _ in trims]           # dt column of each trim switch
    ibias_x = next(slot)                           # ibias bus -> array ibias bus
    pair_x = [next(slot) for _ in segs]            # s column of each pair
    os_x = [next(slot) for _ in oss]               # offset-segment cascodes

    # ---- y plan (going down from the array) ------------------------------
    y_on = geo["bottom"] - 1500                    # outn bus (m4) centre
    y_op = y_on - 1500                             # outp / ibias buses (m3)
    y_ring_top = y_op - BUS_H / 2 - 1500           # ring sits below the buses
    in_top = y_ring_top - 410 - 600                # ring inner edge
    # cascode row: bottoms aligned, so the gate bars share one cbias line
    wf = {n: int(round(dac.cascode_w(k) * 1000))
          for n, k in dac.segments() + dac.os_segments()}
    y_cb = in_top - 800 - M1H_S / 2 - 400 - max(wf.values())
    y_cbias = y_cb - HEAD_GAP - HEAD_H / 2
    # a p-tap strip between the rows: the tallest cascodes are 13 um, and
    # LU.2 wants every diffusion within 15 um of a tap
    y_tap1 = y_cbias - M1H_S / 2 - 700
    y_tap0 = y_tap1 - RING_TAP
    w_max = int(round(max(dac.switch_w(k) for _, k in dac.segments()) * 1000))
    y_dev = y_tap0 - 700 - w_max                   # switch diffusion bottom
    y_rd = in_top - 200 - RD_HH                    # resistor centre
    y_vdd = y_rd - RD_TERM                         # vdd bus (m2) centre
    in_bot = min(y_rd - RD_HH, y_dev - HEAD_GAP - HEAD_H) - 400
    y_pin = in_bot - 600 - 410 - 3000

    # ---- trim switches: [dt | trim_k | ibias] ----------------------------
    for net, x in zip(trims, trim_x):
        x0 = x - col_x(0)
        cols, gates, hy, yc = device_row(d, x0, y_dev, 1)
        up_m2(d, cols[0], yc, bus_y[net] + 200)             # dt: up to its bus
        d.square("via2", cols[0], bus_y[net], VIA2)
        up_m3(d, cols[1], yc, y_op + BUS_H / 2)             # ibias
        gate_down(d, gates[0], hy, y_pin, f"trim{net[2:]}")
    ib_lo = min(x - col_x(0) + col_x(1) for x in trim_x) - M3_W / 2
    ib_hi = max(max(x - col_x(0) + col_x(1) for x in trim_x), ibias_x) + M3_W / 2
    d.rect("m3", ib_lo, y_op - BUS_H / 2, ib_hi, y_op + BUS_H / 2)
    d.rect("m3pin", ib_lo, y_op - BUS_H / 2, ib_lo + 600, y_op + BUS_H / 2)
    d.label("m3lbl", "ibias", ib_lo + 300, y_op)
    # the array's ibias bus: one drop onto this one
    d.square("via2", ibias_x, bus_y["ibias"], VIA2)
    d.rect("m2", ibias_x - M2_W / 2, y_op - 200, ibias_x + M2_W / 2,
           bus_y["ibias"] + 200)
    d.square("via2", ibias_x, y_op, VIA2)

    # ---- per segment: cascode, and under it the switch pair ---------------
    #   drop (m2) -> cascode strap ... cascode drain -> m2 -> pair source
    #   pair: [outp | d | s | db | vdd]   -- unipolar: the other side dumps
    op_xs, dump_xs = [], []
    src = "c" if dac.CASCODE else "s"
    for name, x in zip(segs, pair_x):
        x0 = x - col_x(1)
        sw = int(round(dac.switch_w(dict(dac.segments())[name]) * 1000))
        cols, gates, hy, yc = device_row(d, x0, y_dev, 2, sw)
        up_m3(d, cols[0], yc, y_op + BUS_H / 2)                  # outp
        up_m4(d, cols[2], yc, y_on + BUS_H / 2)                  # dump -> vdd
        if dac.CASCODE:
            y_xv, y_strap = cascode(d, x, y_cb, wf[name])
            up_m2(d, cols[1], yc, y_xv + M2P_H / 2)              # pair source
            d.rect("m2", x - M2_W / 2, y_strap - M2P_H / 2, x + M2_W / 2,
                   bus_y[f"{src}_{name}"] + 200)                 # the drop
        else:
            up_m2(d, cols[1], yc, bus_y[f"{src}_{name}"] + 200)
        d.square("via2", x, bus_y[f"{src}_{name}"], VIA2)
        gate_down(d, gates[0], hy, y_pin, f"d_{name}")
        gate_down(d, gates[1], hy, y_pin, f"db_{name}")
        op_xs.append(cols[0])
        dump_xs.append(cols[2])

    # ---- offset segments: a cascode, and its drain straight up to outp ----
    for name, x in zip(oss, os_x):
        y_xv, y_strap = cascode(d, x, y_cb, wf[name])
        d.rect("m2", x - M2_W / 2, y_strap - M2P_H / 2, x + M2_W / 2,
               bus_y[f"c_{name}"] + 200)                         # the drop
        d.square("via2", x, bus_y[f"c_{name}"], VIA2)
        # drain: m2 down past the cbias line, then up in m3 to the outp bus
        y_turn = y_tap0 - 700
        d.rect("m2", x - M2_W / 2, y_turn - M2P_H / 2, x + M2_W / 2,
               y_xv + M2P_H / 2)
        d.square("via2", x, y_turn, VIA2)
        d.rect("m3", x - M3_W / 2, y_turn - 200, x + M3_W / 2, y_op + BUS_H / 2)
        op_xs.append(x)
    x_last = max(pair_x + os_x)

    ring_x0 = min(vss_x) - 1000                    # ring inner left edge
    if dac.CASCODE:
        # cbias: one m1 line through every gate bar, and an m2 down to its pin
        x_cbp = pair_x[0] - 1500
        d.rect("m1", x_cbp - VIA / 2 - 85, y_cbias - M1H_S / 2,
               x_last + SD + CL + SD / 2, y_cbias + M1H_S / 2)
        d.square("via", x_cbp, y_cbias, VIA)
        d.rect("m2", x_cbp - M2_W / 2, y_pin, x_cbp + M2_W / 2,
               y_cbias + M2P_H / 2)
        d.rect("m2pin", x_cbp - M2_W / 2, y_pin, x_cbp + M2_W / 2, y_pin + 500)
        d.label("m2lbl", "cbias", x_cbp, y_pin + 250)
        # the tap strip, strapped in m1 into the ring's left side
        tx0, tx1 = ring_x0 - 600, x_last + 1500
        d.rect("tap", tx0, y_tap0, tx1, y_tap1)
        d.rect("psdm", tx0 - SDM_ENC, y_tap0 - SDM_ENC, tx1 + SDM_ENC,
               y_tap1 + SDM_ENC)
        d.rect("li", tx0 - RING_TAP / 2, y_tap0 + 40, tx1 - 40, y_tap1 - 40)
        ty = (y_tap0 + y_tap1) / 2
        for x in _array(tx0 + 120, tx1 - 120, LIC, LIC_P):
            d.square("licon", x, ty, LIC)
        for x in _array(tx0 + 120, tx1 - 120, LIC, MCON_P):
            d.square("mcon", x, ty, LIC)
        d.rect("m1", tx0 - RING_TAP / 2, y_tap0, tx1, y_tap1)

    # ---- load resistors, right of everything else --------------------------
    # Flattened into this cell, so the vias below sit on its terminals without
    # stacking contacts across the hierarchy.
    x_rd = max(max(dump_xs), x_last + SD + CL + SD / 2) + 3000 + RD_HW
    inst = top.insert(kdb.CellInstArray(rd.cell_index(),
                                        kdb.Trans(int(x_rd), int(y_rd))))
    inst.flatten()
    for dx in (-RD_X, RD_X):
        for y in (y_rd + RD_TERM, y_rd - RD_TERM):
            d.rect("m1", x_rd + dx - M1P_W / 2, y - M1P_H / 2,
                   x_rd + dx + M1P_W / 2, y + M1P_H / 2)
            d.square("via", x_rd + dx, y, VIA)
    up_m3(d, x_rd - RD_X, y_rd + RD_TERM, y_op + BUS_H / 2)     # R_P -> outp
    op_xs.append(x_rd - RD_X)
    # R_N -> outn: vn's precharge path, and a plain pin -- it stops below the
    # buses, so it never meets the vdd dump bus on m4
    y_outn = y_ring_top + 700
    up_m3(d, x_rd + RD_X, y_rd + RD_TERM, y_outn + 300)
    d.rect("m3pin", x_rd + RD_X - M3_W / 2, y_outn - 300, x_rd + RD_X + M3_W / 2,
           y_outn + 300)
    d.label("m3lbl", "outn", x_rd + RD_X, y_outn)
    # both bottom terminals on the vdd bus, and vdd down to its pin
    x_vr = x_rd + RD_HW + 400                      # vdd riser to the dump bus
    d.rect("m2", x_rd - RD_X - M2P_W / 2, y_vdd - BUS_H / 2,
           x_vr + M2P_W / 2, y_vdd + BUS_H / 2)
    d.rect("m2", x_rd + RD_X - M2_W / 2, y_pin, x_rd + RD_X + M2_W / 2, y_vdd)
    d.rect("m2pin", x_rd + RD_X - M2_W / 2, y_pin, x_rd + RD_X + M2_W / 2,
           y_pin + 500)
    d.label("m2lbl", "vdd", x_rd + RD_X, y_pin + 250)
    d.square("via2", x_vr, y_vdd, VIA2)
    d.rect("m3", x_vr - M3_W / 2, y_vdd - 400, x_vr + M3_W / 2, y_vdd + 400)
    d.square("via3", x_vr, y_vdd, VIA3)
    d.rect("m4", x_vr - M3_W / 2, y_vdd - 200, x_vr + M3_W / 2, y_on + BUS_H / 2)
    dump_xs.append(x_vr)
    # the resistors' own guard ring is li only: strap it to the strip's ring
    gx, gy = x_rd + 1000, y_rd - RD_B
    d.square("mcon", gx, gy, LIC)
    d.rect("m1", gx - 145, in_bot - 600 - 205, gx + 145, gy + 145)

    # ---- output and dump buses -----------------------------------------------
    d.rect("m3", min(op_xs) - M3_W / 2, y_op - BUS_H / 2, max(op_xs) + M3_W / 2,
           y_op + BUS_H / 2)
    d.label("m3lbl", "outp", max(op_xs), y_op)
    d.rect("m3pin", max(op_xs) - M3_W / 2, y_op - BUS_H / 2, max(op_xs) + M3_W / 2,
           y_op + BUS_H / 2)
    d.rect("m4", min(dump_xs) - M3_W / 2, y_on - BUS_H / 2, max(dump_xs) + M3_W / 2,
           y_on + BUS_H / 2)

    # ---- p-tap ring around the strip, and its vss drops -------------------
    inner = kdb.Box(int(min(vss_x) - 1000), int(in_bot),
                    int(x_vr + 800), int(in_top))
    rbox, _ = ring(d, inner, "psdm", "vss")
    for x in vss_x:
        d.square("via2", x, bus_y["vss"], VIA2)
        vy = rbox[3] - RING_TAP / 2                 # ring side centre line
        d.rect("m2", x - M2_W / 2, vy - VIA / 2 - 85, x + M2_W / 2,
               bus_y["vss"] + 200)
        d.square("via", x, vy, VIA)

    full = d.bbox(list(LAYERS))
    d.rect("bound", full.left, full.bottom, full.right, full.top)
    return ly, top


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("-o", "--out", default=os.path.join(here, CELL + ".gds"))
    ap.add_argument("--png", help="also render the layout to this PNG")
    args = ap.parse_args()
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    ly, top = build(out_dir)
    ly.write(args.out)
    bb = top.dbbox()
    print(f"wrote {args.out}: {bb.width():.1f} x {bb.height():.1f} um")
    if args.png:
        render(args.out, args.png)


if __name__ == "__main__":
    main()
