#!/usr/bin/env python3
"""Size the I-DAC unit current source from sky130 Monte-Carlo mismatch.

Next item 4 of silicon/README.md.  The tolerance budget asks for
``dac_sigma_lsb``: the rms error of the 1-LSB element, in LSB.  The behavioural
model draws every element with sigma * sqrt(weight), i.e. independent unit
sources, so the budget is a *relative* unit-current mismatch:

    sigma(I_unit) / I_unit  <=  dac_sigma_lsb            (read from the spec)

For each candidate W/L this finds the gate voltage that gives the unit current
at nominal (tt, no mismatch), then draws ``--draws`` Monte-Carlo samples of
``--per-draw`` identical devices at that bias (tt_mm).  Deviation from the
within-draw mean is pure local mismatch -- tt_mm carries no die-level term (a
pair probe measured a correlation of 0.03) -- so every device is one sample.

    ./osic.ps1 python3 ngspice/idac_unit_size.py          # ~15 min
    ./osic.ps1 python3 ngspice/idac_unit_size.py --quick  # ~5 min

Writes ngspice/idac_unit_size.json next to this file.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import subprocess
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
SPEC = HERE.parent / "spec" / "fdvpd_spec.json"
LIB = "/foss/pdks/sky130A/libs.tech/ngspice/sky130.lib.spice"

#: drain voltage of the unit source.  The cascode above it holds the source's
#: drain roughly constant; 0.5 V is where the stack puts it (see the report).
VDS = 0.5

#: candidate (W, L) in um.  Long and narrow: at 195 nA a short device would sit
#: deep in weak inversion, where gm/Id is highest and V_T mismatch hurts most.
CANDIDATES = [(1.0, 8.0), (1.0, 16.0), (0.5, 8.0), (0.5, 16.0), (1.0, 32.0),
              (0.42, 24.0), (2.0, 32.0)]


def ngspice(netlist: str) -> str:
    with tempfile.NamedTemporaryFile("w", suffix=".spice", delete=False) as fh:
        fh.write(netlist)
        path = fh.name
    out = subprocess.run(["ngspice", "-b", path], capture_output=True,
                         text=True, timeout=1800)
    return out.stdout + out.stderr


def values(text: str, name: str) -> list[float]:
    return [float(v) for v in re.findall(rf"^{name}\s*=\s*(\S+)", text, re.M)]


def bias_for(w: float, l: float, i_unit: float) -> dict:
    """Gate voltage giving i_unit at nominal, plus V_dsat and gm/Id there."""
    text = ngspice(f"""* bias search
.lib {LIB} tt
XM d g 0 0 sky130_fd_pr__nfet_01v8 W={w} L={l}
Vd d 0 {VDS}
Vg g 0 0.6
.control
dc Vg 0.3 1.2 0.002
let id = -i(Vd)
meas dc vg_at when id={i_unit}
let gm = deriv(id)
let gmid = gm / id
meas dc gmid_at find gmid when id={i_unit}
.endc
.end
""")
    vg = values(text, "vg_at")
    gmid = values(text, "gmid_at")
    if not vg:
        raise RuntimeError(f"no bias point for W={w} L={l}:\n{text[-800:]}")
    # V_dsat from the saturation knee: the Vds at which Id reaches 98 % of its
    # value at VDS, with the gate held at the bias just found.
    knee = ngspice(f"""* knee
.lib {LIB} tt
XM d g 0 0 sky130_fd_pr__nfet_01v8 W={w} L={l}
Vd d 0 {VDS}
Vg g 0 {vg[0]}
.control
dc Vd 0 {VDS} 0.002
let id = -i(Vd)
let ref = id[length(id) - 1]
meas dc vdsat when id={0.98 * i_unit}
.endc
.end
""")
    vdsat = values(knee, "vdsat")
    return {"vg": vg[0], "gm_id": gmid[0] if gmid else float("nan"),
            "vdsat_98": vdsat[0] if vdsat else float("nan")}


def mismatch_all(sizes: list[tuple[float, float, float]], draws: int,
                 per_draw: int) -> list[float]:
    """Relative sigma for every (W, L, Vg) in ``sizes``, from one MC loop.

    Every size sits in the same netlist, each with its own gate source, so one
    ``mc_source`` -- which reloads the whole sky130 library, ~10 s -- draws all
    of them at once.  Local mismatch is independent per instance, so sharing a
    draw costs nothing statistically.
    """
    lines, reads, norms = [], [], []
    for s_i, (w, l, vg) in enumerate(sizes):
        lines.append(f"Vg{s_i} g{s_i} 0 {vg}")
        base = s_i * draws * per_draw
        for j in range(per_draw):
            lines.append(f"XM{s_i}_{j} d{s_i}_{j} g{s_i} 0 0 "
                         f"sky130_fd_pr__nfet_01v8 W={w} L={l}")
            lines.append(f"V{s_i}_{j} d{s_i}_{j} 0 {VDS}")
            reads.append(f"  let s[{base}+k*{per_draw}+{j}] = -i(V{s_i}_{j})")
        idx = [f"s[{base}+k*{per_draw}+{j}]" for j in range(per_draw)]
        norms.append(f"  let m = ({' + '.join(idx)}) / {per_draw}")
        norms += [f"  let {x} = {x} / m - 1" for x in idx]
    n = draws * per_draw
    prints = "\n".join(
        f"let v{s_i} = s[{s_i*n},{(s_i+1)*n - 1}]\n"
        f"let rel{s_i} = sqrt(mean(v{s_i} * v{s_i}) * {per_draw} / {per_draw - 1})\n"
        f"print rel{s_i}" for s_i in range(len(sizes)))
    body = "\n".join(lines)
    loop = "\n".join(reads + norms)
    text = ngspice(f"""* unit mismatch, all sizes per draw
.lib {LIB} tt_mm
{body}
.control
let s = unitvec({len(sizes) * n})
let k = 0
while k < {draws}
  mc_source
  op
{loop}
  let k = k + 1
end
{prints}
.endc
.end
""")
    out = []
    for s_i in range(len(sizes)):
        v = values(text, f"rel{s_i}")
        if not v:
            raise RuntimeError(f"MC failed:\n{text[-1500:]}")
        out.append(v[0])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--draws", type=int, default=100)
    ap.add_argument("--per-draw", type=int, default=8)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    if args.quick:
        args.draws = 20

    spec = json.loads(SPEC.read_text())
    i_unit = spec["operating_point"]["idac"]["unit_current_a"]
    budget = spec["tolerance_budget"]["dac_sigma_lsb"]["limit"]
    n = args.draws * args.per_draw
    print(f"unit current {i_unit*1e9:.1f} nA, budget sigma <= {budget:.3g} "
          f"(relative), {n} samples per size\n")
    print(f"{'W':>5} {'L':>5} {'area':>6} {'Vg':>6} {'gm/Id':>6} "
          f"{'Vdsat':>6} {'sigma':>7} {'A*sig':>7} {'margin':>7}")

    biases = [bias_for(w, l, i_unit) for w, l in CANDIDATES]
    sigmas = mismatch_all([(w, l, b["vg"]) for (w, l), b in
                           zip(CANDIDATES, biases)], args.draws, args.per_draw)
    rows = []
    for (w, l), b, sig in zip(CANDIDATES, biases, sigmas):
        # standard error of a sigma estimate from n samples
        se = sig / math.sqrt(2 * (n - 1))
        row = {"w_um": w, "l_um": l, "area_um2": w * l, **b,
               "sigma_rel": sig, "sigma_se": se,
               "pelgrom_um_pct": 100 * sig * math.sqrt(w * l),
               "margin": budget / sig}
        rows.append(row)
        print(f"{w:5.2f} {l:5.1f} {w*l:6.1f} {b['vg']:6.3f} {b['gm_id']:6.1f} "
              f"{b['vdsat_98']:6.3f} {100*sig:6.2f}% {row['pelgrom_um_pct']:6.2f} "
              f"{row['margin']:6.2f}x")

    out = HERE / "idac_unit_size.json"
    out.write_text(json.dumps({"i_unit_a": i_unit, "budget_sigma_rel": budget,
                               "vds_v": VDS, "samples_per_size": n,
                               "rows": rows}, indent=2))
    print(f"\nwrote {out.relative_to(HERE.parent)}")


if __name__ == "__main__":
    main()
