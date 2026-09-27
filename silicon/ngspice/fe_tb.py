#!/usr/bin/env python3
"""FDVPD front end, end to end: DAC + ramp + switches + C_SAR, many cycles.

    .\\osic.ps1 python3 ngspice/fe_tb.py

The architecture check before any SAR or control is designed on top of it.
Each 20 ns reference cycle, with ideal timing:

    REF + 0.5 ns   the DAC takes the next code (vp is disconnected)
    REF + 10 ns    S1 closes: vp settles to the DAC's outp, vn to VDD via R_N
    REF - 0.3 ns   S1 opens: both sides of C_SAR hold
    REF            S2 closes: the ramp sink discharges vn at SR
    REF + dt       S2 opens; vn - vp is the residue, read 1 ns later

The unipolar DAC puts (c + 154) I_u R_D across the pair, and the ramp takes
SR * dt off it, so the residue is zero at dt = c * T_PD / 1024 + 300 ps.  The
codes follow a fractional channel (steps of 32, wrapping) and each dt is the
ideal lock point plus a known offset of up to +-20 ps, so a fit of residue
against code and dt separates gain, offset and linearity, and what is left is
what the previous cycle leaves behind: DAC settling, vn's recovery, charge
injection.

Ideal bias sources for now (ibias 50 uA, cbias 1.30 V, the ramp's two 200 uA);
C_SAR is 0.5 pF of ideal capacitance per side.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import subprocess
import sys
import tempfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
LAYOUT = HERE.parent / "layout"
sys.path.insert(0, str(LAYOUT))
import fe_netlist as fe  # noqa: E402
import idac_netlist as dac  # noqa: E402

LIB = "/foss/pdks/sky130A/libs.tech/ngspice/sky130.lib.spice"
VDD = 1.8
T_REF = 20e-9
T_PD = 2e-9
T_OS = 300e-12
EDGE = 0.1e-9
PFET = "sky130_fd_pr__pfet_01v8"


def code_bits(code):
    msb, lsb = code >> 6, code & 63
    on = {}
    for name, _ in dac.segments():
        on[name] = int(msb > int(name[1:])) if name.startswith("t") else \
            (lsb >> int(name[1:])) & 1
    return on


#: diagnosis switch: FE_IDEAL=s1 / s2 / s1s2 swaps those PMOS switches for ideal
#: ones (200 Ohm on, no charge injection, no body effect)
IDEAL = os.environ.get("FE_IDEAL", "")


def switches():
    out = [".model swi sw vt=0.9 vh=0 ron=200 roff=1e12"]
    if "s1" in IDEAL:
        # ideal: closed when the (active-low) gate is below 0.9 V
        out += ["Ss1p outp vp vs1i 0 swi", "Ss1n outn vn vs1i 0 swi",
                f"Bs1i vs1i 0 v={VDD} - v(s1)"]
    else:
        # the front end's own encode switches; here s1 is the active-low gate
        out += [ln.replace(" s1 ", " s1b ").replace(" s1_n ", " s1 ")
                for ln in fe.encode()]
    if "s2" in IDEAL:
        out += ["Ss2 rs vn vs2i 0 swi", f"Bs2i vs2i 0 v={VDD} - v(s2)",
                "Ss2d rs vdd s2 0 swi"]
    else:
        # the front end's own steering switch and gate bias, on s2_n = s2
        out += [ln.replace("s2_n", "s2") for ln in fe.steer()] + ["Isw vdd vg 20u"]
    return out


def pwl(points):
    """points: [(t, v)] already sorted; returns a PWL string."""
    return " ".join(f"{t:.9e} {v:g}" for t, v in points)


def square(windows, lo, hi, t_end):
    """Level `hi` inside each (t0, t1) window, `lo` outside."""
    pts = [(0.0, lo)]
    for t0, t1 in windows:
        pts += [(t0, lo), (t0 + EDGE, hi), (t1, hi), (t1 + EDGE, lo)]
    pts.append((t_end, lo))
    return pwl(pts)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--start", type=int, default=100, help="first DAC code")
    ap.add_argument("--step", type=int, default=32,
                    help="code step per cycle (121 walks the whole range, as det_cosim)")
    ap.add_argument("--cycles", type=int, default=14)
    ap.add_argument("--tag", default="", help="suffix for the JSON name")
    args = ap.parse_args()
    rng = np.random.default_rng(1)
    n_cyc = args.cycles
    codes = [(args.start + args.step * k) % 1024 for k in range(n_cyc)]
    offsets = rng.uniform(-20e-12, 20e-12, n_cyc)
    offsets[0] = 0.0
    t0 = 25e-9                       # first REF: the first encode is long
    refs = [t0 + k * T_REF for k in range(n_cyc)]
    dts = [c * T_PD / 1024 + T_OS + o for c, o in zip(codes, offsets)]
    t_end = refs[-1] + 5e-9

    lines = [f"* front end, {n_cyc} cycles", f".lib {LIB} tt",
             dac.netlist(), (LAYOUT / "ramp_sink_ref.spice").read_text(),
             ".option rshunt=1e12",
             f"Vdd vdd 0 {VDD}", "Vss vss 0 0",
             "Ib vdd ib 50u", "Vcb cbias 0 1.30",
             "Irr vdd nbias 200u", "Irc vdd rcbias 200u",
             f"Xdut {' '.join(dac.ports()).replace('ibias', 'ib')} idac",
             "Xr vss rs nbias rcbias ramp_sink",
             # S1: encode switches (PMOS, on when low); S2: the ramp switch
             *switches(),
             "Cp vp 0 0.5p", "Cn vn 0 0.5p",
             ".save v(vp) v(vn)"]
    for k in range(len(dac.REF_TRIM)):
        lines.append(f"Vtrim{k} trim{k} 0 {VDD if 64 >> k & 1 else 0}")
    # DAC codes: code k applied at REF_{k-1} + 0.5 ns (the first from t = 0)
    for name, _ in dac.segments():
        pts_d, pts_db = [(0.0, code_bits(codes[0])[name] * VDD)], []
        for k in range(1, n_cyc):
            t = refs[k - 1] + 0.5e-9
            prev, cur = code_bits(codes[k - 1])[name], code_bits(codes[k])[name]
            pts_d += [(t, prev * VDD), (t + EDGE, cur * VDD)]
        pts_d.append((t_end, pts_d[-1][1]))
        pts_db = [(t, VDD - v) for t, v in pts_d]
        lines.append(f"Vd_{name} d_{name} 0 pwl({pwl(pts_d)})")
        lines.append(f"Vdb_{name} db_{name} 0 pwl({pwl(pts_db)})")
    # S1 on (low) from REF_{k-1} + 10 ns (first: from 0) to REF_k - 0.3 ns
    enc = [(0.0, refs[0] - 0.3e-9)] + [(r + 10e-9, r + T_REF - 0.3e-9) for r in refs[:-1]]
    lines.append(f"Vs1 s1 0 pwl({square(enc, VDD, 0.0, t_end)})")
    # S2 (s2_n) low during each ramp: the current steers into vn
    ramps = [(r, r + dt) for r, dt in zip(refs, dts)]
    lines.append(f"Vs2 s2 0 pwl({square(ramps, VDD, 0.0, t_end)})")
    lines += [".control", "set wr_singlescale", "set wr_vecnames",
              f"tran 5p {t_end} 0 20p", "wrdata /tmp/fe.txt v(vp) v(vn)",
              ".endc", ".end"]
    with tempfile.NamedTemporaryFile("w", suffix=".spice", delete=False) as fh:
        fh.write("\n".join(lines) + "\n")
    out = subprocess.run(["ngspice", "-b", fh.name], capture_output=True,
                         text=True, timeout=7200)
    try:
        rows = [ln.split() for ln in open("/tmp/fe.txt") if ln.strip()]
    except OSError:
        raise RuntimeError(out.stdout[-3000:] + out.stderr[-2000:])
    data = np.array([[float(x) for x in r] for r in rows[1:]])
    t, vp, vn = data[:, 0], data[:, 1], data[:, 2]
    res = []
    for r, dt in zip(refs, dts):
        i = np.searchsorted(t, r + dt + 1e-9)
        res.append(float(vn[i] - vp[i]))
    res = np.array(res)
    # fit residue = a*(c + 154) - s*dt + b; skip cycle 0 (no previous cycle)
    A = np.column_stack([np.array(codes[1:]) + dac.N_OS, -np.array(dts[1:]),
                         np.ones(n_cyc - 1)])
    coef, *_ = np.linalg.lstsq(A, res[1:], rcond=None)
    fit = A @ coef
    err = res[1:] - fit
    lsb_v = coef[0]
    print(f"cycles {n_cyc - 1} (first skipped)")
    print(f"DAC step       {lsb_v*1e6:.1f} uV per code (ideal 781.25)")
    print(f"ramp slope     {coef[1]*1e-9:.4f} mV/ps (ideal 0.4000)")
    print(f"gain ratio     {lsb_v / (coef[1] * T_PD / 1024):.4f} "
          f"(DAC step over ramp's per-code time; 1 = on target)")
    print(f"offset         {coef[2]*1e3:.2f} mV")
    print(f"residual rms   {np.std(err)*1e6:.1f} uV = {np.std(err)/coef[1]*1e15:.0f} fs, "
          f"max {np.max(np.abs(err))*1e6:.1f} uV "
          f"(SAR LSB 195 uV; DAC LSB 781 uV)")
    for c, dt, e in zip(codes[1:], dts[1:], err):
        print(f"  code {c:4d}  dt {dt*1e12:7.1f} ps  residual {e*1e6:8.1f} uV")
    tag = "_".join(x for x in (IDEAL, args.tag) if x)
    (HERE / f"fe_tb{'_' + tag if tag else ''}.json").write_text(json.dumps({
        "codes": codes, "dt_s": dts, "residue_v": res.tolist(),
        "fit": {"lsb_v": lsb_v, "slope_v_per_s": coef[1], "offset_v": coef[2]},
        "residual_v": err.tolist()}, indent=2))


if __name__ == "__main__":
    main()
