#!/usr/bin/env python3
"""FDVPD ramp sink -- sky130A layout generator.

Draws the wide-swing cascode ramp sink that ngspice/tb_ramp_casc.spice proved
(branch C), minus the things that are not part of this cell:

  * I_R and I_b are ideal sources in the testbench; here they are the two
    current INPUTS `nbias` and `cbias`.  Where they come from is Next item 2.
  * C_SAR is the SAR's capacitor array, so `outp` is a port and the 1 pF is
    somebody else's layout.

    device  W/L         nf x W_f   S     G      D      notes
    Mr      20 / 1      4 x 5      vss   nbias  nrc    matched pair, ABBA
    M1c     20 / 1      4 x 5      vss   nbias  nxc    matched pair, ABBA
    Mrc     20 / 1      4 x 5      nrc   cbias  nbias  cascode, reference
    M2c     20 / 1      4 x 5      nxc   cbias  outp   cascode, ramp
    Mwb      5 / 1      1 x 5      vss   cbias  cbias  wide-swing bias
    Msw     20 / 0.15   4 x 5      vdd   pre    outp   precharge (pfet)
    Mdum     5 / 1      2 x 5      vss   vss    vss    dummies around Mr/M1c

Mr and M1c are the only pair whose matching the testbench measures (the 3.7 mV
`mirror_match`), so they are the ones drawn common-centroid: one diffusion,
shared sources, drains in A B B A order, with a grounded dummy finger at each
end so every active finger sees the same poly neighbourhood.

Floorplan, bottom to top -- each row has its S/D routing tracks (m2) below it
and its gate bar (m1) above it, so the two never cross:

    [ n-ring / nwell ]  row 3  Msw
    [ p-ring         ]  row 2  Mrc | M2c | Mwb
                        row 1  Mdum Mr M1c M1c Mr Mdum
    routing channel on the right: one m3 vertical per inter-row net

Everything is on a 5 nm grid; the numbers below are in nm.

    python3 ramp_sink.py            # writes ramp_sink.gds next to this file
    python3 ramp_sink.py --png out  # also renders a picture (needs klayout.lay)
"""
import argparse
import os

import klayout.db as kdb

CELL = "ramp_sink"

# sky130 GDS layer / datatype
LAYERS = {
    "nwell": (64, 20), "diff": (65, 20), "tap": (65, 44),
    "poly": (66, 20), "licon": (66, 44), "npc": (95, 20),
    "nsdm": (93, 44), "psdm": (94, 20),
    "li": (67, 20), "mcon": (67, 44),
    "m1": (68, 20), "m1lbl": (68, 5), "m1pin": (68, 16), "via": (68, 44),
    "m2": (69, 20), "m2lbl": (69, 5), "via2": (69, 44),
    "m3": (70, 20), "m3lbl": (70, 5), "m3pin": (70, 16),
    "bound": (235, 4),
}

# --- geometry (nm) -----------------------------------------------------------
SD = 360          # S/D column width: licon 170 + 95 to each gate (licon.11 55)
LIC = 170         # licon / mcon edge
LIC_P = 340       # licon pitch (licon.2 spacing 170)
MCON_P = 360      # mcon pitch (ct.2 spacing 190)
VIA = 150         # via1 edge
VIA2 = 200        # via2 edge
M1W = 320         # m1 strip: via1 150 + 2 x 85 (via.5a)
M2W = 320         # m2 track
M3W = 330         # m3 vertical: via2 200 + 2 x 65
POLY_END = 130    # poly.8 endcap
GB_GAP = 250      # diff top -> gate-bar poly (keeps poly licon 350 off diff)
GB_H = 370        # gate-bar poly: licon 170 + 2 x 100 (npc enclosure)
TRACK_P = 600     # m2 track pitch
TRACK_0 = 1000    # diff bottom -> first track centre (clears dummy heads)
CHAN_P = 800      # m3 pitch in the routing channel
RING_TAP = 410    # tap ring width: licon 170 + 2 x 120 (licon.7)
RING_GAP = 600    # ring inner edge -> anything inside
SDM_ENC = 125     # nsdm / psdm enclosure of diff and tap
NWELL_ENC = 180   # nwell enclosure of pdiff / ntap


def snap(v):
    return int(round(v / 5.0)) * 5


class Drawer:
    def __init__(self, layout, cell):
        self.ly, self.cell = layout, cell
        self.idx = {k: layout.layer(*v) for k, v in LAYERS.items()}

    def rect(self, layer, x0, y0, x1, y1):
        x0, x1 = sorted((snap(x0), snap(x1)))
        y0, y1 = sorted((snap(y0), snap(y1)))
        self.cell.shapes(self.idx[layer]).insert(kdb.Box(x0, y0, x1, y1))
        return (x0, y0, x1, y1)

    def square(self, layer, cx, cy, a):
        return self.rect(layer, cx - a / 2, cy - a / 2, cx + a / 2, cy + a / 2)

    def label(self, layer, text, x, y):
        self.cell.shapes(self.idx[layer]).insert(
            kdb.Text(text, kdb.Trans(snap(x), snap(y))))

    def bbox(self, layers):
        box = kdb.Box()
        for name in layers:
            box += self.cell.bbox_per_layer(self.idx[name])
        return box


def _array(lo, hi, size, pitch):
    """Centres of the most `size` squares at `pitch` that fit in [lo, hi]."""
    n = int((hi - lo - size) // pitch) + 1
    span = (n - 1) * pitch + size
    start = lo + snap((hi - lo - span) / 2) + size / 2
    return [start + i * pitch for i in range(n)]


# =============================================================================
# a row of multi-finger transistors
# =============================================================================
class Row:
    """One row of devices sharing a gate net.

    `diffs` is a list of (cols, dummies): `cols` the S/D net of each column
    (nf + 1 of them), `dummies` the finger indices whose gate goes to vss
    instead of the row's gate net.
    """

    def __init__(self, d, x0, yb, w, l, gate, diffs, gap=800):
        self.d, self.yb, self.w, self.l, self.gate = d, yb, w, l, gate
        self.cols = []        # (x_centre, net, diff_index, col_index)
        self.fingers = []     # (x0, x1, dummy, outer_col_x)
        self.diff_boxes = []
        x = x0
        for di, (cols, dummies) in enumerate(diffs):
            dx0 = x
            for ci, net in enumerate(cols):
                self.cols.append((x + SD / 2, net, di, ci))
                x += SD
                if ci < len(cols) - 1:
                    fi = ci
                    dummy = fi in dummies
                    outer = None
                    if dummy:
                        outer = dx0 + SD / 2 if fi == 0 else x + l + SD / 2
                    self.fingers.append((x, x + l, dummy, outer))
                    x += l
            self.diff_boxes.append((dx0, x))
            x += gap
        self.x1 = x - gap
        self.nets = []
        for _, net, _, _ in self.cols:
            if net not in self.nets:
                self.nets.append(net)
        self.gb0 = yb + w + GB_GAP                  # gate-bar poly bottom
        self.gate_y = self.gb0 + GB_H / 2           # gate-bar centre line

    def track_y(self, net):
        return self.yb - TRACK_0 - self.nets.index(net) * TRACK_P

    @property
    def bottom(self):
        return self.track_y(self.nets[-1]) - M2W / 2

    @property
    def top(self):
        return self.gb0 + GB_H

    def draw(self, sdm, x_right, dummy_net="vss"):
        d, yb, w = self.d, self.yb, self.w
        yt = yb + w
        for dx0, dx1 in self.diff_boxes:
            d.rect("diff", dx0, yb, dx1, yt)
            d.rect(sdm, dx0 - SDM_ENC, yb - SDM_ENC, dx1 + SDM_ENC, yt + SDM_ENC)

        # --- S/D columns: licon -> li -> mcon -> m1 strip down to its track
        lic_y = _array(yb + 60, yt - 60, LIC, LIC_P)
        li_y0, li_y1 = lic_y[0] - LIC / 2 - 80, lic_y[-1] + LIC / 2 + 80
        mc_y = _array(li_y0, li_y1, LIC, MCON_P)
        dummy_outer = {f[3] for f in self.fingers if f[2]}
        head_top = yb - 300
        head_bot = head_top - GB_H
        for cx, net, _, _ in self.cols:
            for y in lic_y:
                d.square("licon", cx, y, LIC)
            li_bot = head_bot + 20 if any(abs(cx - o) < 1 for o in dummy_outer) else li_y0
            d.rect("li", cx - LIC / 2, li_bot, cx + LIC / 2, li_y1)
            for y in mc_y:
                d.square("mcon", cx, y, LIC)
            ty = self.track_y(net)
            d.rect("m1", cx - M1W / 2, ty - VIA / 2 - 85, cx + M1W / 2,
                   mc_y[-1] + LIC / 2 + 60)
            d.square("via", cx, ty, VIA)

        # --- gates
        active = [f for f in self.fingers if not f[2]]
        for fx0, fx1, dummy, outer in self.fingers:
            if not dummy:
                d.rect("poly", fx0, yb - POLY_END, fx1, self.gate_y)
                continue
            # dummy: poly down to a contact head, li across to the outer column
            fc = (fx0 + fx1) / 2
            d.rect("poly", fx0, head_bot, fx1, yt + POLY_END)
            hx0, hx1 = min(fx0, fc - 185), max(fx1, fc + 185)
            d.rect("poly", hx0, head_bot, hx1, head_top)
            ly = head_bot + 100 + LIC / 2
            d.square("licon", fc, ly, LIC)
            d.rect("npc", hx0, head_bot, hx1, head_top)
            d.rect("li", min(outer - LIC / 2, fc - LIC / 2), head_bot + 20,
                   max(outer + LIC / 2, fc + LIC / 2), head_bot + 350)

        # --- gate bar: poly -> licon -> li -> mcon -> m1, out to the channel
        bx0, bx1 = active[0][0], active[-1][1]
        gb0 = self.gb0
        d.rect("poly", bx0, gb0, bx1, gb0 + GB_H)
        d.rect("npc", bx0, gb0, bx1, gb0 + GB_H)
        for x in _array(bx0 + 100, bx1 - 100, LIC, LIC_P):
            d.square("licon", x, self.gate_y, LIC)
        d.rect("li", bx0 + 20, gb0 + 20, bx1 - 20, gb0 + GB_H - 20)
        for x in _array(bx0 + 50, bx1 - 50, LIC, MCON_P):   # m1.4: 30 past bar
            d.square("mcon", x, self.gate_y, LIC)
        d.rect("m1", bx0 + 20, self.gate_y - M1W / 2, x_right, self.gate_y + M1W / 2)

    def tracks(self, x0, x1, ring_x=None):
        """m2 tracks.  vss / vdd run out onto the guard ring and via down to it."""
        for net in self.nets:
            y = self.track_y(net)
            if ring_x and net in ring_x[2]:
                # past the ring's centre line, so the via keeps its enclosure
                e = VIA / 2 + 85
                self.d.rect("m2", ring_x[0] - e, y - M2W / 2, ring_x[1] + e, y + M2W / 2)
                for rx in ring_x[:2]:
                    self.d.square("via", rx, y, VIA)
            else:
                self.d.rect("m2", x0, y - M2W / 2, x1, y + M2W / 2)


def stack(d, x, y, from_m1=True):
    """m1 (optional) -> via1 -> m2 pad -> via2 -> m3, centred on (x, y)."""
    if from_m1:
        d.square("via", x, y, VIA)
    d.square("m2", x, y, 400)
    d.square("via2", x, y, VIA2)


def ring(d, box, tap_sdm, net):
    """Substrate / well tap guard ring just outside `box`, strapped in li + m1."""
    x0, y0 = box.left - RING_GAP - RING_TAP, box.bottom - RING_GAP - RING_TAP
    x1, y1 = box.right + RING_GAP + RING_TAP, box.top + RING_GAP + RING_TAP
    t = RING_TAP
    sides = [(x0, y0, x1, y0 + t), (x0, y1 - t, x1, y1),
             (x0, y0, x0 + t, y1), (x1 - t, y0, x1, y1)]
    for sx0, sy0, sx1, sy1 in sides:
        d.rect("tap", sx0, sy0, sx1, sy1)
        d.rect("li", sx0 + 40, sy0 + 40, sx1 - 40, sy1 - 40)
        d.rect("m1", sx0, sy0, sx1, sy1)
    e = SDM_ENC
    for sx0, sy0, sx1, sy1 in sides:
        d.rect(tap_sdm, sx0 - e, sy0 - e, sx1 + e, sy1 + e)
    # contacts: full length along top and bottom, corners excluded on the sides
    for yc in (y0 + t / 2, y1 - t / 2):
        for x in _array(x0 + 120, x1 - 120, LIC, LIC_P):
            d.square("licon", x, yc, LIC)
        for x in _array(x0 + 120, x1 - 120, LIC, MCON_P):
            d.square("mcon", x, yc, LIC)
    for xc in (x0 + t / 2, x1 - t / 2):
        for y in _array(y0 + t + 200, y1 - t - 200, LIC, LIC_P):
            d.square("licon", xc, y, LIC)
        for y in _array(y0 + t + 200, y1 - t - 200, LIC, MCON_P):
            d.square("mcon", xc, y, LIC)
    d.label("m1lbl", net, x0 + t / 2, y0 + t / 2)
    d.rect("m1pin", x0, y0, x0 + t, y0 + t)
    return (x0, y0, x1, y1), (x0 + t / 2, x1 - t / 2)


def build():
    ly = kdb.Layout()
    ly.dbu = 0.001
    top = ly.create_cell(CELL)
    d = Drawer(ly, top)

    W, L_N, L_P = 5000, 1000, 150

    # ---- row 1: Mdum Mr M1c M1c Mr Mdum, one diffusion, shared sources
    r1_cols = ["vss", "vss", "nrc", "vss", "nxc", "vss",
               "nxc", "vss", "nrc", "vss", "vss"]
    r1 = Row(d, 0, 0, W, L_N, "nbias", [(r1_cols, {0, 9})])

    # ---- row 2: Mrc | M2c | Mwb, three diffusions under one cbias gate bar
    r2_diffs = [(["nrc", "nbias", "nrc", "nbias", "nrc"], set()),
                (["nxc", "outp", "nxc", "outp", "nxc"], set()),
                (["vss", "cbias"], set())]
    probe = Row(d, 0, 0, W, L_N, "cbias", r2_diffs)
    r2_yb = r1.top + 600 + TRACK_0 + (len(probe.nets) - 1) * TRACK_P + M2W / 2
    r2 = Row(d, 0, snap(r2_yb), W, L_N, "cbias", r2_diffs)

    # ---- routing channel, right of both rows
    chan_nets = ["nrc", "nxc", "nbias", "cbias", "outp", "pre"]
    chan0 = max(r1.x1, r2.x1) + 1000
    chan = {n: chan0 + i * CHAN_P for i, n in enumerate(chan_nets)}
    x_end = chan[chan_nets[-1]] + 400

    for r in (r1, r2):
        r.draw("nsdm", chan[r.gate] + 200)
    stack(d, chan["nbias"], r1.gate_y)
    stack(d, chan["cbias"], r2.gate_y)
    # nrc / nxc / nbias / cbias / outp land on their row-2 tracks from m3
    for net in ("nrc", "nxc"):
        stack(d, chan[net], r1.track_y(net), from_m1=False)
        # internal, not ports -- named so the extracted netlist can be probed
        d.label("m2lbl", net, r1.cols[0][0], r1.track_y(net))
    for net in ("nrc", "nbias", "nxc", "outp", "cbias"):
        stack(d, chan[net], r2.track_y(net), from_m1=False)

    # p-ring: sized from everything drawn so far, then vss tracks reach it
    inner = d.bbox(["diff", "poly", "li", "m1", "m2"])
    inner = kdb.Box(inner.left, inner.bottom, max(inner.right, x_end), inner.top)
    pring_box, pring_x = ring(d, inner, "psdm", "vss")
    for r in (r1, r2):
        r.tracks(r.cols[0][0] - 400, x_end, (pring_x[0], pring_x[1], {"vss"}))

    # ---- row 3: precharge pfet in its own nwell + n-ring
    r3_cols = ["vdd", "outp", "vdd", "outp", "vdd"]
    probe = Row(d, 0, 0, W, L_P, "pre", [(r3_cols, set())])
    r3_ring_bot = pring_box[3] + 1000 + NWELL_ENC + SDM_ENC
    r3_yb = (r3_ring_bot + RING_TAP + RING_GAP
             + TRACK_0 + (len(probe.nets) - 1) * TRACK_P + M2W / 2)
    r3 = Row(d, 0, snap(r3_yb), W, L_P, "pre", [(r3_cols, set())])
    r3.draw("psdm", chan["pre"] + 200)
    stack(d, chan["pre"], r3.gate_y)
    stack(d, chan["outp"], r3.track_y("outp"), from_m1=False)

    n_inner = kdb.Box(inner.left, snap(r3.bottom), inner.right, snap(r3.top))
    nring_box, nring_x = ring(d, n_inner, "nsdm", "vdd")
    r3.tracks(r3.cols[0][0] - 400, x_end, (nring_x[0], nring_x[1], {"vdd"}))
    d.rect("nwell", nring_box[0] - NWELL_ENC, nring_box[1] - NWELL_ENC,
           nring_box[2] + NWELL_ENC, nring_box[3] + NWELL_ENC)

    # ---- m3 verticals; the four signal ports leave through the top edge
    y_top = nring_box[3] + NWELL_ENC + 600
    spans = {
        "nrc": (r1.track_y("nrc"), r2.track_y("nrc")),
        "nxc": (r1.track_y("nxc"), r2.track_y("nxc")),
        "nbias": (r1.gate_y, y_top),
        "cbias": (r2.track_y("cbias"), y_top),
        "outp": (r2.track_y("outp"), y_top),
        "pre": (r3.gate_y, y_top),
    }
    for net, (ya, yb) in spans.items():
        x = chan[net]
        y0, y1 = min(ya, yb) - 200, max(ya, yb) + (0 if yb == y_top else 200)
        d.rect("m3", x - M3W / 2, y0, x + M3W / 2, y1)
        if yb == y_top:
            d.rect("m3pin", x - M3W / 2, y_top - 500, x + M3W / 2, y_top)
            d.label("m3lbl", net, x, y_top - 250)

    full = d.bbox(list(LAYERS))
    d.rect("bound", full.left, full.bottom, full.right, full.top)
    return ly, top


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("-o", "--out", default=os.path.join(here, CELL + ".gds"))
    ap.add_argument("--png", help="also render the layout to this PNG")
    args = ap.parse_args()
    ly, top = build()
    ly.write(args.out)
    bb = top.dbbox()
    print(f"wrote {args.out}: {bb.width():.2f} x {bb.height():.2f} um")
    if args.png:
        render(args.out, args.png)


def render(gds, png):
    import klayout.lay as klay
    view = klay.LayoutView()
    view.load_layout(gds)
    view.max_hier()
    view.set_config("background-color", "#ffffff")
    view.set_config("grid-visible", "false")
    colours = {
        "nwell": 0xa0a0ff, "diff": 0x40c040, "tap": 0x208020, "poly": 0xff4040,
        "li": 0xa050ff, "m1": 0x4080ff, "m2": 0xff9f40, "m3": 0x40c0c0,
    }
    it = view.begin_layers()
    while not it.at_end():
        lp = it.current()
        name = next((k for k, v in LAYERS.items()
                     if v == (lp.source_layer, lp.source_datatype)), None)
        props = lp.dup()
        props.visible = name in colours
        if name in colours:
            props.fill_color = props.frame_color = colours[name]
            props.dither_pattern = 5 if name in ("m2", "m3", "nwell") else 1
            props.transparent = True
        view.set_layer_properties(it, props)
        it.next()
    view.zoom_fit()
    view.save_image(png, 1400, 1800)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
