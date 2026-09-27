#!/usr/bin/env python3
"""FDVPD top level -- cell ``fdvpd``: the detector's four analog cells and the
control macro, placed and routed.

    [ idac ....................................... ]
    ======  gap: DAC switch / trim strips (m2) straight down to the macro,
            outn / outp / cbias / ibias leave on m4 tracks
    [ fdvpd_ctrl (LibreLane) ]  | channel |  [ sarfe  (mirrored) ]  | vss (m4)
                                |  m2     |  [ ramp_sink         ]  | vdd (m5)
                                |  cols   |  [ bias              ]  |  trunks

The three analog cells are placed with the M135 mirror, (x, y) -> (-y, -x),
so the pins on their top edges face west into the channel and their power
rings -- all on their left sides -- face up into the gap above each cell.
The channel is two-layer: every net owns one m2 column and runs m3 along its
pins' own y, so no two nets can cross on one layer.  The four nets that start
in the DAC take m4 instead, top track to the easternmost column, so their
drops never cross each other's tracks.

Power: via stacks on each ring in the gaps, vss on m4 lines to an m4 trunk,
vdd on m5 lines to an m5 trunk; the macro's own m5 straps and the DAC's dump
bus (vdd) and ring (vss) join the same two trunks.

Reads the cells from build/<cell>/<cell>.gds (verify_cell.sh writes them) and
the macro from the newest LibreLane run, so run those first:

    python3 top.py -o build/fdvpd/fdvpd.gds --png build/fdvpd/fdvpd.png
"""
import argparse
import glob
import os
import sys

import klayout.db as kdb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "digital"))
import pins_def  # noqa: E402
from ramp_sink import LAYERS, Drawer  # noqa: E402

LAYERS.update({"capm": (89, 44), "m2lbl": (69, 5), "m2pin": (69, 16), "via3": (70, 44),
               "m4": (71, 20), "m4lbl": (71, 5), "m4pin": (71, 16),
               "via4": (71, 44), "m5": (72, 20), "m5lbl": (72, 5),
               "m5pin": (72, 16)})

CELL = "fdvpd"
MACRO = "fdvpd_ctrl"

# --- widths (nm) --------------------------------------------------------------
W2, W3, W4 = 280, 300, 400        # signal wires on m2 / m3 / m4
COL_P = 800                        # m2 column pitch in the channel
COL4_P = 1000                      # m4 column pitch
W4_PWR, W5_PWR = 600, 1600         # power lines
TRUNK4_W, TRUNK5_W = 2000, 3200

# --- floorplan (nm) -----------------------------------------------------------
DAC_BOT = -57280                   # idac's bottom edge
MACRO_TOP = -63000                 # 5.7 um gap for the strips and m4 tracks
Y_JOG = -58000                     # m2 jog of the DAC strips, just under the DAC
#: m4 tracks in the DAC-macro gap, top first; source x ordering is the same
GAP_TRACKS = {"outn": -58300, "outp": -59300, "cbias": -60300, "ibias": -61300}
CH_X0 = 206000                     # first channel column
CELL_GAP = 4700                    # between stacked cells: their power lines
#: macro east pin -> sarfe pin
CTRL = (["s1_n", "s2_n", "cmp_clk", "cmp_p", "cmp_n"]
        + [f"bpos{k}" for k in range(6)] + [f"bneg{k}" for k in range(6)]
        + [f"bmid{k}" for k in range(6)])
TAPS = ["vrp20", "vrp10", "vrp5", "vrm", "vrn5", "vrn10", "vrn20"]
#: bias pin -> (cell, pin) for the channel nets between the analog cells
LOCAL = ([("rs", ("ramp_sink", "outp"), ("sarfe", "rs")),
          ("vg", ("bias", "isw"), ("sarfe", "vg")),
          ("nbias", ("bias", "iramp"), ("ramp_sink", "nbias")),
          ("rcbias", ("bias", "icbr"), ("ramp_sink", "cbias"))]
         + [(t, ("bias", t), ("sarfe", t)) for t in TAPS])


def newest_run():
    runs = sorted(glob.glob(os.path.join(HERE, "..", "digital", "runs", "RUN_*")))
    for r in reversed(runs):
        g = os.path.join(r, "final", "gds", MACRO + ".gds")
        if os.path.exists(g) and os.path.exists(os.path.join(r, "final", "pnl", MACRO + ".pnl.v")):
            return r
    raise FileNotFoundError("no LibreLane run with final/gds")


def load(ly, path, name):
    src = kdb.Layout()
    src.read(path)
    c = ly.create_cell(name)
    c.copy_tree(src.cell(name))
    labels = {}
    for li in src.layer_indexes():
        info = src.get_info(li)
        for s in src.cell(name).shapes(li).each():
            if s.is_text():
                labels.setdefault(s.text_string, []).append(
                    (info.layer, s.text_pos.x, s.text_pos.y))
    return c, src, labels


def ring_box(src, name, x, y):
    """Outer bbox of the merged m1 shape under the label at (x, y)."""
    li = src.find_layer(*LAYERS["m1"])
    m1 = kdb.Region(src.cell(name).begin_shapes_rec(li)).merged()
    hit = m1.interacting(kdb.Region(kdb.Box(x - 5, y - 5, x + 5, y + 5)))
    assert not hit.is_empty(), f"{name}: no m1 under ({x}, {y})"
    return hit.bbox()


class Top:
    def __init__(self):
        self.ly = kdb.Layout()
        self.ly.dbu = 0.001
        self.top = self.ly.create_cell(CELL)
        self.d = Drawer(self.ly, self.top)

    # ---- primitives ---------------------------------------------------------
    # wires carry square end caps, half a width past each end point, so two
    # meeting at a corner leave no notch (met4.1 flagged every bare L)
    def hwire(self, layer, x0, x1, y, w):
        self.d.rect(layer, min(x0, x1) - w / 2, y - w / 2, max(x0, x1) + w / 2,
                    y + w / 2)

    def vwire(self, layer, x, y0, y1, w):
        self.d.rect(layer, x - w / 2, min(y0, y1) - w / 2, x + w / 2,
                    max(y0, y1) + w / 2)

    def via2(self, x, y):              # m2 <-> m3
        self.d.square("m2", x, y, 370)
        self.d.square("via2", x, y, 200)
        self.d.square("m3", x, y, 500)

    def via3(self, x, y):              # m3 <-> m4
        self.d.square("m3", x, y, 500)
        self.d.square("via3", x, y, 200)
        self.d.square("m4", x, y, 500)

    def stack(self, x, y, top):
        """Power stack m1 -> m4 (top=4) or m5 (top=5): 2 x 2 vias to m4."""
        d = self.d
        d.square("m1", x, y, 800)
        for dx in (-160, 160):
            for dy in (-160, 160):
                d.square("via", x + dx, y + dy, 150)
        d.square("m2", x, y, 800)
        for dx in (-200, 200):
            for dy in (-200, 200):
                d.square("via2", x + dx, y + dy, 200)
                d.square("via3", x + dx, y + dy, 200)
        d.square("m3", x, y, 800)
        d.square("m4", x, y, 1200)
        if top == 5:
            d.square("via4", x, y, 800)
            d.square("m5", x, y, 1600)

    def label(self, layer, text, x, y, w, h):
        self.d.rect(layer + "pin", x - w / 2, y - h / 2, x + w / 2, y + h / 2)
        self.d.label(layer + "lbl", text, x, y)


def build():
    t = Top()
    ly, top, d = t.ly, t.top, t.d

    # ---- cells ------------------------------------------------------------------
    build_dir = os.path.join(HERE, "build")
    run = newest_run()
    cells = {}
    for name in ("idac", "sarfe", "ramp_sink", "bias"):
        cells[name] = load(ly, os.path.join(build_dir, name, name + ".gds"), name)
    cells[MACRO] = load(ly, os.path.join(run, "final", "gds", MACRO + ".gds"), MACRO)

    M135 = 7                          # kdb fixpoint code: (x, y) -> (-y, -x)
    # the mirrored cell's west edge sits at CELL_X, stacked from the top down
    n_m2 = len(CTRL) + len(LOCAL)
    CELL_X = CH_X0 + (n_m2 - 1) * COL_P + 2000 + 3 * COL4_P + 2000
    place = {"idac": kdb.Trans(0, 0), MACRO: kdb.Trans(0, MACRO_TOP - 110000)}
    y_top = DAC_BOT - CELL_GAP         # sarfe's top edge
    gaps = {}
    for name in ("sarfe", "ramp_sink", "bias"):
        bb = cells[name][1].cell(name).bbox()
        # mirrored bbox: x in [-top, -bottom], y in [-right, -left]
        dx = CELL_X + bb.top
        dy = y_top + bb.left
        place[name] = kdb.Trans(M135, False, dx, dy)
        gaps[name] = (y_top, y_top + CELL_GAP)          # the gap above it
        y_top = dy - bb.right - CELL_GAP
    for name, (c, _, _) in cells.items():
        top.insert(kdb.CellInstArray(c.cell_index(), place[name]))

    def pin(cell, net, layer=None):
        """Placed position of a cell's top-level label."""
        hits = [(l, x, y) for l, x, y in cells[cell][2].get(net, [])
                if layer is None or l == layer]
        assert hits, f"{cell}: no label {net}"
        p = place[cell] * kdb.Point(hits[0][1], hits[0][2])
        return p.x, p.y

    # ---- DAC switch and trim strips: macro north pins up to the DAC ----------
    macro_pins = pins_def.rows()       # (name, layer, x, y, side), macro coords
    mx, my = 0, MACRO_TOP - 110000
    for name, layer, x, y, side in macro_pins:
        if side != "N":
            continue
        base, k = name.split("[")[0], int(name.split("[")[1][:-1])
        dac_net = f"trim{k}" if base == "trim" else \
            ("d_" if base == "dsw" else "db_") + pins_def.SEG[k]
        xd, _ = pin("idac", dac_net)
        xm = mx + x
        t.vwire("m2", xm, my + y, Y_JOG, W2)
        t.hwire("m2", xm, xd, Y_JOG, W2)
        t.vwire("m2", xd, Y_JOG, DAC_BOT + 250, W2)

    # ---- the four DAC nets: m4 tracks in the gap, then m4 columns -------------
    # the top track's net takes the easternmost column
    col4 = {"outn": CELL_X - 2000, "outp": CELL_X - 2000 - COL4_P,
            "cbias": CELL_X - 2000 - 2 * COL4_P,
            "ibias": CELL_X - 2000 - 3 * COL4_P}
    for net, yt in GAP_TRACKS.items():
        xs, ys = pin("idac", net)
        if net == "cbias":                               # an m2 pin: drop, then up
            t.vwire("m2", xs, DAC_BOT + 250, yt, 320)
            t.via2(xs, yt)
            t.via3(xs, yt)
        else:                                            # an m3 pin inside the DAC
            t.via3(xs, ys)
            t.vwire("m4", xs, ys, yt, W4)
        t.hwire("m4", xs, col4[net], yt, W4)

    # ---- channel: each net one m2 column (or m4 for the DAC nets) -------------
    def to_cell(x_col, cell, net):
        """m3 from a column to a mirrored cell's west-edge pin."""
        px, py = pin(cell, net, 70)
        t.hwire("m3", x_col, px, py, W3)
        return py

    col = CH_X0
    east = {n: (x, y) for n, l, x, y, s in macro_pins if s == "E"}
    for net in CTRL:
        pname = net if not net[-1].isdigit() else f"{net[:-1]}[{net[-1]}]"
        ex, ey = east[pname]
        ym = my + ey
        t.hwire("m3", mx + ex, col, ym, W3)
        t.via2(col, ym)
        ys = to_cell(col, "sarfe", net)
        t.vwire("m2", col, ym, ys, W2)
        t.via2(col, ys)
        col += COL_P
    for net, a, b in LOCAL:
        ya = to_cell(col, *a)
        yb = to_cell(col, *b)
        t.via2(col, ya)
        t.via2(col, yb)
        t.vwire("m2", col, ya, yb, W2)
        col += COL_P
    assert col - COL_P < col4["ibias"] - 1500, "channel columns overlap the m4 ones"
    targets = {"outn": ("sarfe", "outn"), "outp": ("sarfe", "outp"),
               "cbias": ("bias", "dcbias"), "ibias": ("bias", "idac")}
    for net, (cell, p) in targets.items():
        yt = GAP_TRACKS[net]
        px, py = pin(cell, p, 70)
        t.vwire("m4", col4[net], yt, py, W4)
        t.via3(col4[net], py)
        t.hwire("m3", col4[net], px, py, W3)

    # ---- power ------------------------------------------------------------------
    right = max((place[n] * cells[n][1].cell(n).bbox()).right
                for n in ("sarfe", "ramp_sink", "bias"))
    x_ss = right + 2000 + TRUNK4_W // 2
    x_dd = x_ss + TRUNK4_W // 2 + 2000 + TRUNK5_W // 2
    ss_ys, dd_ys = [], []

    def ring_stacks(cell, net, y_line, top_layer):
        """Stubs from the cell's ring (its mirrored top side) up to y_line."""
        _, src, labels = cells[cell]
        _, lx, ly_ = [v for v in labels[net] if v[0] == 68][0]
        b = ring_box(src, cell, lx, ly_)
        # the ring's left side, mirrored: y = -b.left, spanning x = -b.top .. -b.bottom
        p0 = place[cell] * kdb.Point(b.left, b.top)
        p1 = place[cell] * kdb.Point(b.left, b.bottom)
        x0, x1 = sorted((p0.x, p1.x))
        y_ring = p0.y                  # the ring's outer edge
        span = x1 - x0 - 4000
        n = max(2, span // 20000 + 1)
        xs = [x0 + 2000 + i * span // (n - 1) for i in range(n)]
        for x in xs:
            t.vwire("m1", x, y_ring, y_line, 600)     # the cap reaches 0.3 into the ring
            t.stack(x, y_line, top_layer)
        return min(xs)

    for cell in ("sarfe", "ramp_sink", "bias"):
        g0, g1 = gaps[cell]
        y_ss, y_dd = g1 - 1300, g1 - 3300
        x = ring_stacks(cell, "vss", y_ss, 4)
        t.hwire("m4", x, x_ss, y_ss, W4_PWR)
        ss_ys.append(y_ss)
        if "vdd" in cells[cell][2] and any(v[0] == 68 for v in cells[cell][2]["vdd"]):
            x = ring_stacks(cell, "vdd", y_dd, 5)
            t.hwire("m5", x, x_dd, y_dd, W5_PWR)
            dd_ys.append(y_dd)

    # DAC vdd: its m4 dump bus, extended east past the vdd riser, up to m5
    y_on, x_riser = -26320, 194980
    t.hwire("m4", x_riser, 197600, y_on, 400)
    t.d.square("m4", 197000, y_on, 1200)
    t.d.square("via4", 197000, y_on, 800)
    t.hwire("m5", 196200, x_dd, y_on, W5_PWR)
    dd_ys.append(y_on)
    # DAC vss: a stub off its strip ring's right side
    x_ring, y_vss = 196790, -42000
    t.hwire("m1", x_ring, 199200, y_vss, 600)
    t.stack(198600, y_vss, 4)
    t.hwire("m4", 198600, x_ss, y_vss, W4_PWR)
    ss_ys.append(y_vss)
    # the macro's m5 straps (VPWR 26.73 - 28.33, VGND 30.03 - 31.63), east
    y_vpwr, y_vgnd = my + 27530, my + 30830
    t.hwire("m5", mx + 198000, x_dd, y_vpwr, W5_PWR)
    # m5 has to enclose the via4 by 0.31 um
    t.hwire("m5", mx + 198000, x_ss, y_vgnd, W5_PWR)
    t.d.square("via4", x_ss, y_vgnd, 800)
    dd_ys.append(y_vpwr)
    ss_ys.append(y_vgnd)
    # trunks
    t.vwire("m4", x_ss, max(ss_ys) + 800, min(ss_ys) - 800, TRUNK4_W)
    t.vwire("m5", x_dd, max(dd_ys) + 800, min(dd_ys) - 800, TRUNK5_W)
    t.label("m4", "vss", x_ss, min(ss_ys), TRUNK4_W, 1000)
    t.label("m5", "vdd", x_dd, min(dd_ys), TRUNK5_W, 1600)

    # ---- the detector's own pins: the macro's south edge -------------------------
    for name, layer, x, y, side in macro_pins:
        if side == "S":
            t.label("m2", name, mx + x, my + y, W2, 1000)

    return ly, top


def render(gds, png):
    """Like ramp_sink.render, metals only and up to m5 (what the top adds)."""
    import klayout.lay as klay
    view = klay.LayoutView()
    view.load_layout(gds)
    view.max_hier()
    view.set_config("background-color", "#ffffff")
    view.set_config("grid-visible", "false")
    colours = {"m1": 0x4080ff, "m2": 0xff9f40, "m3": 0x40c0c0,
               "m4": 0xb040c0, "m5": 0x808080, "capm": 0x40a040}
    it = view.begin_layers()
    while not it.at_end():
        lp = it.current()
        name = next((k for k, v in LAYERS.items()
                     if v == (lp.source_layer, lp.source_datatype)), None)
        props = lp.dup()
        props.visible = name in colours
        if name in colours:
            props.fill_color = props.frame_color = colours[name]
            props.dither_pattern = 5
            props.transparent = True
        view.set_layer_properties(it, props)
        it.next()
    view.zoom_fit()
    view.save_image(png, 1600, 1600)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-o", "--out", default=os.path.join(HERE, "build", CELL, CELL + ".gds"))
    ap.add_argument("--png", help="also render the layout to this PNG")
    args = ap.parse_args()
    ly, top = build()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    ly.write(args.out)
    bb = top.dbbox()
    print(f"wrote {args.out}: {bb.width():.1f} x {bb.height():.1f} um")
    if args.png:
        render(args.out, args.png)


if __name__ == "__main__":
    main()
