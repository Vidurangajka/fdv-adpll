#!/usr/bin/env python3
"""Fit a front-end residue run (fe_tb*.json): linear, and + quadratic terms.

    python3 ngspice/fe_fit.py fe_tb.json [fe_tb_s1.json ...]
"""
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
for name in sys.argv[1:] or ["fe_tb.json"]:
    d = json.loads((HERE / name).read_text())
    c = (np.array(d["codes"][1:]) + 154) / 1000.0
    t = np.array(d["dt_s"][1:]) * 1e9
    r = np.array(d["residue_v"][1:])
    out = [name]
    for label, cols in (("lin", [c, -t, np.ones_like(t)]),
                        ("+dt^2", [c, -t, np.ones_like(t), t ** 2])):
        A = np.column_stack(cols)
        coef, *_ = np.linalg.lstsq(A, r, rcond=None)
        err = r - A @ coef
        s = coef[1] * 1e9
        out.append(f"{label}: rms {err.std()*1e6:6.1f} uV ({err.std()/s*1e15:4.0f} fs)")
        if label == "lin":
            out.append(f"gain {coef[0]*1e-3/(s*2e-9/1024):.4f}")
    print("   ".join(out))
