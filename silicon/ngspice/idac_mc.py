#!/usr/bin/env python3
"""Current DAC Monte-Carlo: segment weights per mismatch draw (sky130 tt_mm).

    .\\osic.ps1 python3 ngspice/idac_mc.py --draws 200      # ~1 h
    python ngspice/idac_mc_loop.py                          # then: through the loop

The DAC's error is the unit sources' mismatch, so each draw is one operating
point of the unit array alone: the reference diode (64 fixed + the 64-unit
trim bit, i.e. the nominal N = 128) fed with I_bias, and every output unit
with its drain held at the cascode-source voltage.  One voltage source per
segment reads that segment's current.  The cascodes and switches are left
out -- the cascode holds the unit drains, so their mismatch only moves Vds by
millivolts, against a 0.9 %/V output conductance.

Every unit is its own instance.  sky130's mismatch model does NOT scale an
``m=`` multiplier like independent devices: an m=64 instance measured the
same relative sigma as one unit (0.77 % vs 0.81 %), i.e. 64 perfectly
correlated copies.  A Monte-Carlo on the m= netlist would overstate segment
mismatch 8x.

Writes ngspice/idac_mc.json: per draw, the 21 segment currents and the
reference diode's gate voltage.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "layout"))
import idac_netlist as dac  # noqa: E402

LIB = "/foss/pdks/sky130A/libs.tech/ngspice/sky130.lib.spice"
NFET = dac.NFET
UNIT = f"W={dac.UNIT_W} L={dac.UNIT_L}"
I_BIAS = 50e-6          # 128 reference units x 390.6 nA
V_CASC_SRC = 0.543      # cascode source at mid-scale, tt 27 C (idac_tb headroom)


def netlist(draws: int, seed: int) -> str:
    L = ["* idac Monte-Carlo, every unit its own instance",
         f".lib {LIB} tt_mm",
         f"Ib 0 ib {I_BIAS}",            # I_bias into the diode node
         ".option rshunt=1e12"]
    n_ref = dac.REF_FIXED + dac.REF_TRIM[-1]
    for j in range(n_ref):
        L.append(f"XD{j} ib ib 0 0 {NFET} {UNIT}")
    segs = dac.segments()
    for k, (name, n) in enumerate(segs):
        L.append(f"Vc{k} c{k} 0 {V_CASC_SRC}")
        for j in range(n):
            L.append(f"XU{k}_{j} c{k} ib 0 0 {NFET} {UNIT}")
    prints = " ".join(f"i(vc{k})" for k in range(len(segs)))
    L += [".control",
          f"setseed {seed}",
          "let k = 0",
          f"while k < {draws}",
          "  mc_source",
          "  op",
          f"  print v(ib) {prints}",
          "  let k = k + 1",
          "end",
          ".endc", ".end"]
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--batch", type=int, default=10,
                    help="draws per ngspice process: mc_source keeps every "
                         "reloaded circuit (~6 MB/s here), and one 200-draw "
                         "process hung the Docker engine")
    args = ap.parse_args()
    segs = dac.segments()
    vg, cur = [], [[] for _ in segs]
    for b in range(0, args.draws, args.batch):
        nb = min(args.batch, args.draws - b)
        with tempfile.NamedTemporaryFile("w", suffix=".spice", delete=False) as fh:
            fh.write(netlist(nb, seed=1 + b))     # a distinct seed per batch
        out = subprocess.run(["ngspice", "-b", fh.name], capture_output=True,
                             text=True, timeout=3600).stdout
        v = [float(x) for x in re.findall(r"^v\(ib\)\s*=\s*(\S+)", out, re.M)]
        c = [[-float(x) for x in re.findall(rf"^i\(vc{k}\)\s*=\s*(\S+)", out, re.M)]
             for k in range(len(segs))]
        if len(v) != nb or any(len(x) != nb for x in c):
            raise RuntimeError(f"batch at {b}: expected {nb}, got {len(v)}\n{out[-2000:]}")
        vg += v
        for k in range(len(segs)):
            cur[k] += c[k]
        print(f"  {len(vg)}/{args.draws} draws", flush=True)
    n = len(vg)
    draws = [{"v_ib": vg[d],
              "segments": {name: cur[k][d] for k, (name, _) in enumerate(segs)}}
             for d in range(n)]
    (HERE / "idac_mc.json").write_text(json.dumps(
        {"i_bias": I_BIAS, "n_ref": dac.REF_FIXED + dac.REF_TRIM[-1],
         "v_casc_src": V_CASC_SRC, "units": dict(segs), "draws": draws}, indent=1))
    # quick summary: relative spread of each segment against its ideal share
    import statistics as st
    for name, units in segs:
        rel = []
        for d in draws:
            tot = sum(d["segments"].values())
            rel.append(d["segments"][name] / (tot * units / 1023) - 1)
        print(f"  {name:4s} {units:3d} units   sigma {st.pstdev(rel)*100:7.4f} %   "
              f"(unit-independent expectation {0.75/units**0.5:.4f} %)")
    print(f"wrote idac_mc.json: {n} draws")


if __name__ == "__main__":
    main()
