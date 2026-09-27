#!/usr/bin/env python3
"""FDVPD SAR front end -- cell ``sarfe``, sky130A layout generator.

Draws what fe_netlist.py describes: the encode and ramp switches, the two CDAC
halves with their three-state bottom-plate switches, and the StrongARM
comparator.  The same approach as ramp_sink.py, generalised to rows whose
fingers each have their own gate net:

  * every device row has its S/D nets on m2 tracks BELOW it (an m1 column
    from each S/D contact down to its net's track) and its gate nets on m2
    tracks ABOVE it (an m1 stub from each finger's head up to its track)
  * a net that leaves its row does so along its track to the ROUTING CHANNEL
    on the right: one m3 vertical per net, which is also where every pin
    is (top edge).  vdd and vss tracks end on their guard ring instead.

    [ n-ring / nwell ]  P    S1p S1n | M5/M6, the reset PFETs, the output
                             inverters' PFETs, the ramp drivers' PFETs
    [ p-ring         ]  C    comparator NMOS: tail, inputs, latch, inverters
                        S    the ramp's steering pair and its drivers' NMOS,
                             the gate-bias replica and diode, vg's MOS decap
                        Nn   n-side bottom-plate switches, weight 1 .. 32
                        Np   p-side bottom-plate switches
                        CDAC n: fixed MIM | unit rows, weight 32 at the bottom
                        CDAC p                                       channel ->

The CDAC's top plates (vp, vn) are the MIM top -- capm, via3 up to an m4
sheet over the whole half -- so the high-impedance side sees the least
parasitic; the bottom plates (m3) are what the switches drive.  Each weight's
units share one m3 plate, whose lead runs right to a via2 and on in m2 to
its channel line.

    python3 sarfe.py -o build/sarfe/sarfe.gds --png build/sarfe/sarfe.png
"""
import argparse
import os
import sys

import klayout.db as kdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fe_netlist as fe  # noqa: E402
from ramp_sink import (LAYERS, Drawer, render, ring, snap,  # noqa: E402
                       LIC, LIC_P, MCON_P, VIA, VIA2, SDM_ENC, RING_TAP,
                       NWELL_ENC, _array)

LAYERS.update({"lvtn": (125, 44), "m2lbl": (69, 5), "m2pin": (69, 16), "via3": (70, 44),
               "m4": (71, 20), "m4lbl": (71, 5), "m4pin": (71, 16),
               "capm": (89, 44)})

CELL = "sarfe"

# --- device geometry (nm) ------------------------------------------------------
L = 150                    # every device in the cell is minimum length
SD = 360                   # S/D column: licon 170 + 95 to each gate
HEAD_W, HEAD_H = 290, 370  # poly gate head: licon 170 + 60 / + 100
HEAD_GAP = 250             # tallest diffusion -> heads
M1P_W = 260                # m1 on a gate head / stub: via1 150 + 2 x 55
M1C_W = 320                # m1 S/D column
M2_W = 320                 # m2 track: via1 150 + 2 x 85
M3_W = 330                 # m3 / m4 line: via 200 + 2 x 65
TRACK_P = 600
SD_TRACK_0 = 800           # diffusion bottom -> first S/D track centre
G_TRACK_0 = 600            # head pad top -> first gate track centre
DIFF_GAP = 800             # between diffusions in a row
ROW_GAP = 1200             # between one row's top track and the next's bottom
CHAN_P = 800               # m3 pitch in the routing channel

# --- CDAC (nm) -------------------------------------------------------------------
CAP = int(round(fe.UNIT_W * 1000))         # unit capm edge
UNIT_P = CAP + 1000                        # along a weight's row (capm.2a 840)
ROW_P = CAP + 1500                         # capm.11: 1.34 to the next plate
M3_ENC = 140                               # capm.3: m3 around capm
FIX = int(round(fe.FIX_SIDE * 1000))
V3_OFF = 400                               # via3 pair, +- from a unit's centre


class FRow:
    """A row of devices on one baseline, each finger with its own gate net.

    ``diffs``: list of (w, cols, gates[, l]) -- w in nm, cols the S/D net of
    each column (len(gates) + 1 of them), gates the gate net of each finger,
    l the gate length (nm, default L).
    """

    def __init__(self, d, x0, yb, diffs, sdm):
        self.d, self.yb, self.sdm = d, yb, sdm
        self.cols, self.fings, self.boxes = [], [], []
        x = x0
        for w, cols, gates, *ln in diffs:
            gl = ln[0] if ln else L
            assert len(cols) == len(gates) + 1
            dx0 = x
            for i, net in enumerate(cols):
                self.cols.append((x + SD / 2, net, w))
                x += SD
                if i < len(gates):
                    self.fings.append((x + gl / 2, gates[i], gl))
                    x += gl
            self.boxes.append((dx0, x, w))
            x += DIFF_GAP
        self.x1 = x - DIFF_GAP
        self.x0 = x0
        wmax = max(spec[0] for spec in diffs)
        self.head_b = yb + wmax + HEAD_GAP
        self.head_y = self.head_b + HEAD_H / 2
        self.sd_nets = list(dict.fromkeys(n for _, n, _ in self.cols))
        self.g_nets = list(dict.fromkeys(f[1] for f in self.fings))

    def sd_y(self, net):
        return self.yb - SD_TRACK_0 - self.sd_nets.index(net) * TRACK_P

    def g_y(self, net):
        return (self.head_y + 200 + G_TRACK_0
                + self.g_nets.index(net) * TRACK_P)

    @property
    def bottom(self):
        return self.sd_y(self.sd_nets[-1]) - M2_W / 2

    @property
    def top(self):
        return self.g_y(self.g_nets[-1]) + M2_W / 2

    def extent(self):
        """Everything the row draws below m2: diffusion to the m1 ends."""
        return kdb.Box(int(self.x0 - SDM_ENC), int(self.sd_y(self.sd_nets[-1]) - 160),
                       int(self.x1 + SDM_ENC), int(self.g_y(self.g_nets[-1]) + 160))

    def terminals(self, net):
        """x of every S/D column and finger on ``net``, and their track ys."""
        pts = [(x, self.sd_y(net)) for x, n, _ in self.cols if n == net]
        pts += [(x, self.g_y(net)) for x, n, _ in self.fings if n == net]
        return pts

    def draw(self):
        d, yb = self.d, self.yb
        for dx0, dx1, w in self.boxes:
            d.rect("diff", dx0, yb, dx1, yb + w)
            d.rect(self.sdm, dx0 - SDM_ENC, yb - SDM_ENC, dx1 + SDM_ENC,
                   yb + w + SDM_ENC)
        # S/D: licon -> li -> mcon -> m1 column down to the net's track
        for cx, net, w in self.cols:
            lic = _array(yb + 60, yb + w - 60, LIC, LIC_P)
            li0, li1 = lic[0] - LIC / 2 - 80, lic[-1] + LIC / 2 + 80
            for y in lic:
                d.square("licon", cx, y, LIC)
            d.rect("li", cx - LIC / 2, li0, cx + LIC / 2, li1)
            mc = _array(li0, li1, LIC, MCON_P)
            for y in mc:
                d.square("mcon", cx, y, LIC)
            ty = self.sd_y(net)
            d.rect("m1", cx - M1C_W / 2, ty - VIA / 2 - 85, cx + M1C_W / 2,
                   mc[-1] + LIC / 2 + 60)
            d.square("via", cx, ty, VIA)
        # gates: poly up to a head, then licon -> li -> mcon -> m1 stub up to
        # the gate net's track
        hb, hy = self.head_b, self.head_y
        for gx, net, gl in self.fings:
            d.rect("poly", gx - gl / 2, yb - 130, gx + gl / 2, hb + HEAD_H)
            hw = max(gl, HEAD_W)
            d.rect("poly", gx - hw / 2, hb, gx + hw / 2, hb + HEAD_H)
            d.square("licon", gx, hy, LIC)
            d.rect("li", gx - LIC / 2, hy - LIC / 2 - 80, gx + LIC / 2,
                   hy + LIC / 2 + 80)
            d.square("mcon", gx, hy, LIC)
            ty = self.g_y(net)
            d.rect("m1", gx - M1P_W / 2, hy - 200, gx + M1P_W / 2,
                   ty + VIA / 2 + 85)
            d.square("via", gx, ty, VIA)
        f0, f1 = self.fings[0], self.fings[-1]
        d.rect("npc", f0[0] - max(f0[2], HEAD_W) / 2 - 40, hb,
               f1[0] + max(f1[2], HEAD_W) / 2 + 40, hb + HEAD_H)

    def tracks(self, chan, ring_x=None, rail=None):
        """m2 tracks: to the channel, or for ``rail`` left onto the ring."""
        for net in dict.fromkeys(self.sd_nets + self.g_nets):
            ys = ([self.sd_y(net)] if net in self.sd_nets else []) + \
                 ([self.g_y(net)] if net in self.g_nets else [])
            for y in ys:
                xs = [x for x, yy in self.terminals(net) if yy == y]
                if net == rail:
                    e = VIA / 2 + 85
                    self.d.rect("m2", ring_x - e, y - M2_W / 2, max(xs) + 200,
                                y + M2_W / 2)
                    self.d.square("via", ring_x, y, VIA)
                elif net in chan:
                    self.d.rect("m2", min(xs) - 200, y - M2_W / 2,
                                chan[net] + 200, y + M2_W / 2)
                    self.d.square("via2", chan[net], y, VIA2)
                elif len(ys) == 2:
                    # local, but both S/D and gate (a diode): the two tracks
                    # meet on an m3 jumper just right of the net's terminals
                    xj = max(x for x, _ in self.terminals(net)) + 500
                    self.d.rect("m2", min(xs) - 200, y - M2_W / 2, xj + 200,
                                y + M2_W / 2)
                    self.d.square("via2", xj, y, VIA2)
                    if y == ys[0]:
                        self.d.rect("m3", xj - M3_W / 2, ys[0] - 200,
                                    xj + M3_W / 2, ys[1] + 200)
                else:                                   # local to the row
                    self.d.rect("m2", min(xs) - 200, y - M2_W / 2,
                                max(xs) + 200, y + M2_W / 2)


def alternate(a, b, n):
    """n + 1 S/D columns alternating a, b, a, ..."""
    return [a if i % 2 == 0 else b for i in range(n + 1)]


def weight_rows():
    """CDAC rows, bottom to top: (k, units), weight 32 first."""
    return [(k, fe.WEIGHTS[k][0]) for k in range(len(fe.WEIGHTS) - 1, -1, -1)]


def cdac(d, x0, y0, side):
    """One CDAC half; returns its geometry.

    Fixed MIM at x0 (bottom plate vss, lead left to x0 - 1000), unit rows to
    its right.  Every capm gets via3 up to the half's m4 sheet (the top
    plate).  Each weight's m3 plate leads right to a via2 at ``x_v``.
    """
    xu = x0 + FIX + 1500                    # first unit column (capm.11)
    n_max = max(n for _, n in weight_rows())
    x_v = xu + (n_max - 1) * UNIT_P + CAP + 1500
    leads = {}
    for r, (k, n) in enumerate(weight_rows()):
        y = y0 + r * ROW_P
        for i in range(n):
            cx = xu + i * UNIT_P
            d.rect("capm", cx, y, cx + CAP, y + CAP)
            for dx in (-V3_OFF, V3_OFF):
                for dy in (-V3_OFF, V3_OFF):
                    d.square("via3", cx + CAP / 2 + dx, y + CAP / 2 + dy, 200)
        d.rect("m3", xu - M3_ENC, y - M3_ENC,
               xu + (n - 1) * UNIT_P + CAP + M3_ENC, y + CAP + M3_ENC)
        ly = y + CAP / 2
        d.rect("m3", xu, ly - 300, x_v + 200, ly + 300)          # the lead
        d.square("via2", x_v, ly, VIA2)
        d.rect("m2", x_v - 200, ly - M2_W / 2, x_v + 200, ly + M2_W / 2)
        leads[f"b{side}{k}"] = (x_v, ly)
    # the fixed part: one square capm, bottom plate to vss on the left
    fy = y0
    d.rect("capm", x0, fy, x0 + FIX, fy + FIX)
    for x in _array(x0 + 300, x0 + FIX - 300, 200, 800):
        for y in _array(fy + 300, fy + FIX - 300, 200, 800):
            d.square("via3", x, y, 200)
    d.rect("m3", x0 - M3_ENC, fy - M3_ENC, x0 + FIX + M3_ENC, fy + FIX + M3_ENC)
    vx, vy = x0 - 1000, fy + FIX / 2
    d.rect("m3", vx - 200, vy - 300, x0, vy + 300)
    d.square("via2", vx, vy, VIA2)
    # the top plate: one m4 sheet over every capm of the half
    rows_top = y0 + (len(fe.WEIGHTS) - 1) * ROW_P + CAP
    top = max(rows_top, fy + FIX)
    sheet = (x0, y0, xu + (n_max - 1) * UNIT_P + CAP, top)
    d.rect("m4", *sheet)
    return {"leads": leads, "vss": (vx, vy), "sheet": sheet, "x_v": x_v,
            "top": top + M3_ENC, "bottom": y0 - M3_ENC}


def build():
    ly = kdb.Layout()
    ly.dbu = 0.001
    top = ly.create_cell(CELL)
    d = Drawer(ly, top)

    # ---- device rows ----------------------------------------------------------
    def plate_row(side):
        diffs = []
        for k, (n, th, tl) in enumerate(fe.WEIGHTS):
            w = 840 if n <= 2 else 1680
            bp = f"b{side}{k}"
            gh, gl = ((f"bpos{k}", f"bneg{k}") if side == "p"
                      else (f"bneg{k}", f"bpos{k}"))
            diffs.append((w, [th, bp, tl], [gh, gl]))
            diffs.append((w, [bp, "vrm"], [f"bmid{k}"]))
        return diffs

    comp = [(4000, ["t", "vss"], ["cmp_clk"]),
            (8000, ["x1", "t", "x1"], ["vn", "vn"]),
            (8000, ["x2", "t", "x2"], ["vp", "vp"]),
            (2000, ["op", "x1"], ["om"]),
            (2000, ["om", "x2"], ["op"]),
            (500, ["cmp_p", "vss"], ["op"]),
            (500, ["cmp_n", "vss"], ["om"])]
    # the ramp switch (fe_netlist.steer): pair, drivers, replica, decap
    sw_w = int(fe.SW_W * 1000)
    steer = [(sw_w, ["rs", "vn", "rs"], ["g", "g"]),
             (sw_w, ["rs", "vdd", "rs"], ["gb", "gb"]),
             (500, ["g", "vss", "gb"], ["s2_n", "g"]),
             (int(fe.REP_W * 1000), ["vg", "m"], ["vg"]),
             (int(fe.DIODE_W * 1000), alternate("m", "vss", fe.DIODE_M),
              ["m"] * fe.DIODE_M, int(fe.DIODE_L * 1000)),
             (int(fe.CAP_W * 1000), alternate("vss", "vss", fe.CAP_M),
              ["vg"] * fe.CAP_M, int(fe.CAP_L * 1000))]
    # the encode switch's NMOS half (vp only) and its gate inverter: a row of their
    # own, because on the steering row they put its middle past LU.2's 15 um
    enc = [(int(fe.S1_WN * 1000), ["outp", "vp"], ["s1"]),
           (1000, ["s1", "vss"], ["s1_n"])]
    pmos = [(int(fe.S1_W * 1000), ["outp", "vp"], ["s1_n"]),
            (int(fe.S1_W * 1000), ["outn", "vn"], ["s1_n"]),
            (2000, ["s1", "vdd"], ["s1_n"]),
            (2000, ["g", "vg", "gb"], ["s2_n", "g"]),
            (2000, ["op", "vdd", "om"], ["om", "op"]),
            (1000, ["op", "vdd", "om"], ["cmp_clk", "cmp_clk"]),
            (1000, ["x1", "vdd", "x2"], ["cmp_clk", "cmp_clk"]),
            (1000, ["cmp_p", "vdd", "cmp_n"], ["op", "om"])]

    # stack bottom-up: each row's lowest track clears the previous row's top
    def place(diffs, sdm, below):
        probe = FRow(d, 0, 0, diffs, sdm)
        yb = below + ROW_GAP - probe.bottom
        return FRow(d, 0, snap(yb), diffs, sdm)

    r_np = place(plate_row("p"), "nsdm", 0)
    r_nn = place(plate_row("n"), "nsdm", r_np.top)
    r_s = place(steer, "nsdm", r_nn.top)
    r_c = place(comp, "nsdm", r_s.top)
    r_e = place(enc, "nsdm", r_c.top)
    n_rows = [r_np, r_nn, r_s, r_c, r_e]
    for r in n_rows:
        r.draw()
    # the encode NMOS is low-Vt: lvtn over its diffusion, 0.18 around the gate
    ex0, ex1, ew = r_e.boxes[0]
    d.rect("lvtn", ex0 - 180, r_e.yb - 180, ex1 + 180, r_e.yb + ew + 180)

    # p-ring around the NMOS rows (devices and their m1; the m2 tracks cross it)
    inner = kdb.Box()
    for r in n_rows:
        inner += r.extent()
    pbox, px = ring(d, inner, "psdm", "vss")

    # P row above, in its own n-ring and nwell, 1.5 um clear of the p-ring
    probe = FRow(d, 0, 0, pmos, "psdm")
    want = pbox[3] + 1500 + NWELL_ENC + RING_TAP + 600      # its extent's bottom
    r_p = FRow(d, 0, snap(want - probe.extent().bottom), pmos, "psdm")
    r_p.draw()
    p_inner = r_p.extent()
    p_inner = kdb.Box(inner.left, p_inner.bottom, max(inner.right, p_inner.right),
                      p_inner.top)
    nbox, nx = ring(d, p_inner, "nsdm", "vdd")
    d.rect("nwell", nbox[0] - NWELL_ENC, nbox[1] - NWELL_ENC,
           nbox[2] + NWELL_ENC, nbox[3] + NWELL_ENC)

    # ---- CDAC halves below the p-ring -------------------------------------------
    half_h = (len(fe.WEIGHTS) - 1) * ROW_P + CAP + 2 * M3_ENC
    y_n = pbox[1] - 2500 - half_h
    y_p = y_n - 3000 - half_h
    cd_p = cdac(d, inner.left + 1500, y_p, "p")
    cd_n = cdac(d, inner.left + 1500, y_n, "n")

    # ---- the channel ---------------------------------------------------------------
    ports = [p for p in fe.ports() if p not in ("vdd", "vss")]
    rows = n_rows + [r_p]
    used = {}
    for r in rows:
        for net in r.sd_nets + r.g_nets:
            used.setdefault(net, set()).add(id(r))
    internal = [n for n in used if n not in ports and n not in ("vdd", "vss")
                and (len(used[n]) > 1 or n in cd_p["leads"] or n in cd_n["leads"])]
    # the dump switch takes vdd inside the p-ring: a channel line for it,
    # strapped to the n-ring
    chan_nets = ports + internal + ["vdd"]
    x_c0 = max(pbox[2], nbox[2], cd_p["x_v"]) + 2000
    chan = {n: x_c0 + i * CHAN_P for i, n in enumerate(chan_nets)}

    for r in n_rows:
        r.tracks(chan, px[0], "vss")
    r_p.tracks(chan, nx[0], "vdd")

    ys = {n: [] for n in chan_nets}
    y_vs = nbox[1] + RING_TAP / 2                 # the n-ring's bottom side
    d.rect("m2", nx[1] - 160, y_vs - M2_W / 2, chan["vdd"] + 200, y_vs + M2_W / 2)
    d.square("via", nx[1], y_vs, VIA)
    d.square("via2", chan["vdd"], y_vs, VIA2)
    ys["vdd"].append(y_vs)
    for r in rows:
        for net in r.sd_nets:
            if net in chan and not (net == "vdd" and r is r_p):
                ys[net].append(r.sd_y(net))
        for net in r.g_nets:
            if net in chan:
                ys[net].append(r.g_y(net))
    # bottom-plate leads: m2 from the half's via2 to the channel
    for cd in (cd_p, cd_n):
        for net, (x, y) in cd["leads"].items():
            d.rect("m2", x - 200, y - M2_W / 2, chan[net] + 200, y + M2_W / 2)
            d.square("via2", chan[net], y, VIA2)
            ys[net].append(y)
    # top plates: the m4 sheet out to the channel, via3 onto vp / vn
    for cd, net in ((cd_p, "vp"), (cd_n, "vn")):
        sx0, sy0, sx1, sy1 = cd["sheet"]
        y = sy1 - 400
        d.rect("m4", sx1 - 400, y - 200, chan[net] + 200, y + 200)
        d.square("via3", chan[net], y, 200)
        ys[net].append(y)
    # the fixed parts' bottom plates: m2 up to the p-ring's bottom side
    vx = cd_p["vss"][0]
    d.rect("m2", vx - 200, cd_p["vss"][1] - 200, vx + 200,
           pbox[1] + RING_TAP / 2 + 160)
    d.square("via", vx, pbox[1] + RING_TAP / 2, VIA)

    y_top = nbox[3] + NWELL_ENC + 1500
    for net in chan_nets:
        x = chan[net]
        lo = min(ys[net]) - 300
        hi = y_top if net in ports else max(ys[net]) + 300
        d.rect("m3", x - M3_W / 2, lo, x + M3_W / 2, hi)
        if net in ports:
            d.rect("m3pin", x - M3_W / 2, y_top - 600, x + M3_W / 2, y_top)
            d.label("m3lbl", net, x, y_top - 300)

    full = d.bbox(list(LAYERS))
    d.rect("bound", full.left, full.bottom, full.right, full.top)
    return ly, top


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("-o", "--out", default=os.path.join(here, CELL + ".gds"))
    ap.add_argument("--png", help="also render the layout to this PNG")
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    ly, top = build()
    ly.write(args.out)
    bb = top.dbbox()
    print(f"wrote {args.out}: {bb.width():.1f} x {bb.height():.1f} um")
    if args.png:
        render(args.out, args.png)


if __name__ == "__main__":
    main()
