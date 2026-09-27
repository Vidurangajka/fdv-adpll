#!/usr/bin/env python3
"""FDVPD bias block -- cell ``bias``, sky130A layout generator.

Draws layout/bias_ref.spice with the sizes ngspice/bias_tb.py settled on:

    [ n-ring / nwell ]  P3   Mo2  (icbr, 8 fingers) | Mo6 (isw)
                        P2   Mo1  (iramp, 8 fingers)
                        P1   M3 / M4 | Mspu | Mo4 | Mo5 | Mo3
    [ p-ring         ]  N2   Msi | Mcr | Mcb (6 fingers) | Mlb
                        N1   M1 / Mspd | M2 (4 fingers)
                        ladder (16 x res_generic_po)   Rs (2 x high_po)   channel ->

Rows are sarfe.FRow -- S/D tracks below, gate tracks above, shared nets out
to an m3 channel on the right that also carries the pins.  No row is wider
than ~23 um, so nothing is more than 15 um from its ring (LU.2).

The two resistors come from magic's own generators (as idac.py's R_D does)
and are flattened in with their labels removed.  The ladder's 16 units are
strapped in series in m1, alternately along the top and the bottom, which
puts every tap -- vrp20 .. vrn20 -- on the top edge; each climbs in m3 to its
own m2 run to the channel.  Rs is two 8 um halves in parallel: a lone 0.35 um
high-sheet resistor is not rpm-legal as generated.

    python3 bias.py -o build/bias/bias.gds --png build/bias/bias.png
    python3 bias.py --ref build/bias/bias_ref_sized.spice   # LVS reference
"""
import argparse
import os
import subprocess
import sys

import klayout.db as kdb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "ngspice"))
import bias_tb  # noqa: E402
from ramp_sink import (LAYERS, Drawer, render, ring, snap,  # noqa: E402
                       LIC, VIA, VIA2, RING_TAP, NWELL_ENC)
from sarfe import FRow, M2_W, M3_W, CHAN_P, ROW_GAP  # noqa: E402

LAYERS.update({"m2lbl": (69, 5), "m2pin": (69, 16), "via3": (70, 44)})

CELL = "bias"
SIZES = bias_tb.SIZES
RLAD = float(SIZES["RLAD"])
CBL = int(round(float(SIZES["CBL"]) * 1000))
CBM = int(SIZES["CBM"])
PORTS = ["vdd", "vss", "iramp", "icbr", "idac", "dcbias", "vrp20", "vrp10",
         "vrp5", "vrm", "vrn5", "vrn10", "vrn20", "isw"]
#: ladder node between unit k-1 and unit k (k = 1..15), and the two ends
LADDER = ["vrp20", "l1", "l2", "l3", "vrp10", "l5", "vrp5", "l7", "vrm", "l9",
          "vrn5", "l11", "vrn10", "l13", "l14", "l15", "vrn20"]


def alternate(a, b, n):
    """n + 1 S/D columns alternating a, b, a, ..."""
    return [a if i % 2 == 0 else b for i in range(n + 1)]


def magic_res(out_dir, name, device, w, l, nx):
    """A resistor array from magic's generator, as a GDS in out_dir."""
    tcl = os.path.join(out_dir, f"{name}.tcl")
    with open(tcl, "w") as fh:
        fh.write(f"""load {name}
box 0 0 0 0
set p [dict merge [sky130::{device}_defaults] {{w {w} l {l} nx {nx}}}]
sky130::{device}_draw $p
gds write {name}.gds
quit -noprompt
""")
    rc = os.path.join(os.environ["PDK_ROOT"], "sky130A", "libs.tech", "magic",
                      "sky130A.magicrc")
    subprocess.run(["magic", "-dnull", "-noconsole", "-rcfile", rc,
                    os.path.basename(tcl)], cwd=out_dir, check=True,
                   capture_output=True)
    return os.path.join(out_dir, f"{name}.gds")


def place_res(ly, top, gds, name, x, y):
    """Flatten a generated resistor into ``top`` at (x, y); return its labels
    (moved) and bbox.  The labels are deleted: two generated cells both say
    R1_0 and B, and magic would join nets by name."""
    src = kdb.Layout()
    src.read(gds)
    sc = src.cell(name)
    lbl = {}
    for li in src.layer_indexes():
        for sh in list(sc.shapes(li).each()):
            if sh.is_text():
                lbl[sh.text_string] = (sh.text_pos.x + x, sh.text_pos.y + y)
                sc.shapes(li).erase(sh)
    c = ly.create_cell(name)
    c.copy_tree(sc)
    inst = top.insert(kdb.CellInstArray(c.cell_index(), kdb.Trans(int(x), int(y))))
    bb = sc.bbox().moved(int(x), int(y))
    inst.flatten()
    return lbl, bb


def pad(d, x, y, layers=("m1",), via=None):
    for layer in layers:
        d.square(layer, x, y, 320)
    if via:
        d.square(via, x, y, VIA if via == "via" else VIA2)


def build(out_dir):
    ly = kdb.Layout()
    ly.dbu = 0.001
    top = ly.create_cell(CELL)
    d = Drawer(ly, top)

    L2, L05, L1 = 2000, 500, 1000
    n1 = [(4000, ["n1", "vss", "st"], ["n1", "n1"], L2),
          (4000, alternate("ns", "pg", 4), ["n1"] * 4, L2)]
    n2 = [(1000, ["pg", "vss"], ["st"], L05),
          (5400, ["dcbias", "ncb"], ["dcbias"], L05),
          (5000, alternate("ncb", "vss", CBM), ["ncb"] * CBM, CBL),
          (4000, ["vrn20", "vss"], ["vrn20"], L1)]
    p1 = [(8000, ["n1", "vdd", "pg"], ["pg", "pg"], L2),
          (420, ["st", "vdd"], ["vss"], 3600),
          (1600, ["dcbias", "vdd"], ["pg"], L2),
          (3200, ["vrp20", "vdd"], ["pg"], L2),
          (8000, ["idac", "vdd", "idac"], ["pg", "pg"], L2)]
    p2 = [(8000, alternate("iramp", "vdd", 8), ["pg"] * 8, L2)]
    p3 = [(8000, alternate("icbr", "vdd", 8), ["pg"] * 8, L2),
          (3200, ["isw", "vdd", "isw"], ["pg", "pg"], L2)]

    def place(diffs, sdm, below):
        probe = FRow(d, 0, 0, diffs, sdm)
        return FRow(d, 0, snap(below + ROW_GAP - probe.extent().bottom), diffs, sdm)

    r_n1 = place(n1, "nsdm", 0)
    r_n2 = place(n2, "nsdm", r_n1.extent().top)
    n_rows = [r_n1, r_n2]
    inner = kdb.Box()
    for r in n_rows:
        r.draw()
        inner += r.extent()
    pbox, px = ring(d, inner, "psdm", "vss")

    below = pbox[3] + 1500 + NWELL_ENC + RING_TAP + 600 - ROW_GAP
    r_p1 = place(p1, "psdm", below)
    r_p2 = place(p2, "psdm", r_p1.extent().top)
    r_p3 = place(p3, "psdm", r_p2.extent().top)
    p_rows = [r_p1, r_p2, r_p3]
    p_inner = kdb.Box()
    for r in p_rows:
        r.draw()
        p_inner += r.extent()
    p_inner = kdb.Box(min(p_inner.left, inner.left), p_inner.bottom,
                      max(p_inner.right, inner.right), p_inner.top)
    nbox, nx = ring(d, p_inner, "nsdm", "vdd")
    d.rect("nwell", nbox[0] - NWELL_ENC, nbox[1] - NWELL_ENC,
           nbox[2] + NWELL_ENC, nbox[3] + NWELL_ENC)

    # ---- resistors, below the p-ring ------------------------------------------------
    lad_gds = magic_res(out_dir, "bias_ladder", "sky130_fd_pr__res_generic_po",
                        1.0, RLAD, 16)
    rs_gds = magic_res(out_dir, "bias_rs", "sky130_fd_pr__res_high_po_0p35",
                       0.35, 8.0, 2)
    # geometry first, from an unplaced copy
    tmp = kdb.Layout()
    tmp.read(lad_gds)
    lad_bb = tmp.cell("bias_ladder").bbox()
    tmp = kdb.Layout()
    tmp.read(rs_gds)
    rs_bb = tmp.cell("bias_rs").bbox()
    y_taps = 1200 + 7 * 700                      # tap runs above the ladder
    y_res = pbox[1] - 2000 - y_taps - max(lad_bb.top, rs_bb.top)
    x_lad = inner.left - lad_bb.left
    lad, lad_box = place_res(ly, top, lad_gds, "bias_ladder", x_lad, y_res)
    x_rs = lad_box.right + 3000 - rs_bb.left
    rs, rs_box = place_res(ly, top, rs_gds, "bias_rs", x_rs, y_res)

    # ---- the channel ------------------------------------------------------------------
    rows = n_rows + p_rows
    used = {}
    for r in rows:
        for net in r.sd_nets + r.g_nets:
            used.setdefault(net, set()).add(id(r))
    ports = [p for p in PORTS if p not in ("vdd",)]
    internal = [n for n in used if n not in PORTS and n != "vdd"
                and (len(used[n]) > 1 or n == "ns")]
    chan_nets = ports + internal
    x_c0 = max(pbox[2], nbox[2], rs_box.right) + 2000
    chan = {n: x_c0 + i * CHAN_P for i, n in enumerate(chan_nets)}
    ys = {n: [] for n in chan_nets}

    # rows: vss / vdd S/D tracks onto their rings, everything else to the
    # channel (a gate on vss -- Mspu's -- takes the channel's vss line)
    for r in n_rows:
        r.tracks(chan, px[0], "vss")
    for r in p_rows:
        r.tracks(chan, nx[0], "vdd")
    for r in rows:
        for net in r.sd_nets:
            if net in chan and not (net == "vss" and r in n_rows):
                ys[net].append(r.sd_y(net))
        for net in r.g_nets:
            if net in chan:
                ys[net].append(r.g_y(net))
    # the channel's vss line joins the p-ring at its bottom-right corner
    y_vr = pbox[1] + RING_TAP / 2               # the ring's bottom side: no tracks there
    d.rect("m2", px[1] - 160, y_vr - M2_W / 2, chan["vss"] + 200, y_vr + M2_W / 2)
    d.square("via", px[1], y_vr, VIA)
    d.square("via2", chan["vss"], y_vr, VIA2)
    ys["vss"].append(y_vr)

    def to_chan(net, x, y, from_m1=True):
        """via1 (from m1) at (x, y), m2 right to the channel line."""
        if from_m1:
            pad(d, x, y, ("m1", "m2"), "via")
        d.rect("m2", x - 160, y - M2_W / 2, chan[net] + 200, y + M2_W / 2)
        d.square("via2", chan[net], y, VIA2)
        ys[net].append(y)

    # ladder: series straps, top for even k, bottom for odd
    for k in range(1, 16):
        side = "R1" if k % 2 == 0 else "R2"
        xa, ya = lad[f"{side}_{k - 1}"]
        xb, _ = lad[f"{side}_{k}"]
        d.rect("m1", xa - 160, ya - 160, xb + 160, ya + 160)
    taps = {"vrp20": lad["R1_0"], "vrn20": lad["R1_15"]}
    for k in (4, 6, 8, 10, 12):
        xa, ya = lad[f"R1_{k - 1}"]
        xb, _ = lad[f"R1_{k}"]
        taps[LADDER[k]] = ((xa + xb) / 2, ya)
    y_run0 = lad_box.top + 1200
    for i, net in enumerate(["vrp20", "vrp10", "vrp5", "vrm", "vrn5", "vrn10", "vrn20"]):
        x, y = taps[net]
        y_run = y_run0 + i * 700
        pad(d, x, y, ("m1", "m2"), "via")
        d.square("m2", x, y, 400)                    # via2 enclosure, stacked
        d.square("via2", x, y, VIA2)
        d.rect("m3", x - M3_W / 2, y - 200, x + M3_W / 2, y_run + 200)
        d.square("via2", x, y_run, VIA2)
        d.rect("m2", x - 200, y_run - M2_W / 2, chan[net] + 200, y_run + M2_W / 2)
        d.square("via2", chan[net], y_run, VIA2)
        ys[net].append(y_run)
    # Rs: the halves in parallel -- top to ns, bottom to vss
    for side, net in (("R1", "ns"), ("R2", "vss")):
        (xa, ya), (xb, _) = rs[f"{side}_0"], rs[f"{side}_1"]
        d.rect("m1", xa - 160, ya - 160, xb + 160, ya + 160)
        to_chan(net, xb, ya)
    # both generated guard rings are li: mcon up at their B label, to vss
    for lbl in (lad, rs):
        x, y = lbl["B"]
        d.square("mcon", x, y, LIC)
        to_chan("vss", x, y)

    y_top = nbox[3] + NWELL_ENC + 1500
    for net in chan_nets:
        x = chan[net]
        lo = min(ys[net]) - 300
        hi = y_top if net in PORTS else max(ys[net]) + 300
        d.rect("m3", x - M3_W / 2, lo, x + M3_W / 2, hi)
        if net in PORTS:
            d.rect("m3pin", x - M3_W / 2, y_top - 600, x + M3_W / 2, y_top)
            d.label("m3lbl", net, x, y_top - 300)

    full = d.bbox(list(LAYERS))
    d.rect("bound", full.left, full.bottom, full.right, full.top)
    return ly, top


def lvs_ref(path):
    """bias_ref.spice with the bench's sizes substituted."""
    with open(path, "w") as fh:
        fh.write(bias_tb.netlist(SIZES))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-o", "--out", default=os.path.join(HERE, CELL + ".gds"))
    ap.add_argument("--png", help="also render the layout to this PNG")
    ap.add_argument("--ref", help="write the sized LVS reference here and stop")
    args = ap.parse_args()
    if args.ref:
        lvs_ref(args.ref)
        return
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
