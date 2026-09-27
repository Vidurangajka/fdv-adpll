#!/usr/bin/env python3
"""Characterise the current DAC netlist (layout/idac_netlist.py) in ngspice.

    .\\osic.ps1 python3 ngspice/idac_tb.py              # tt, 27 C
    .\\osic.ps1 python3 ngspice/idac_tb.py --corner ss --temp 85

Three measurements, each against a number the spec or the model fixes:

  static    every code 0..1023; LSB, full scale, INL / DNL (endpoint fit),
            output common mode, and the lowest output voltage -- the ramp
            window below it starts at 1.4 V
  settling  full-scale step into C_SAR = 1 pF per side; time to 0.1 % and the
            equivalent single-pole tau, against dac_settle_tau <= 4 ns
  R_D       measured, since the poly heads make it non-linear in L

The devices here are nominal -- per-segment ``m=`` multiples -- so INL/DNL is
the systematic part only (finite output impedance, switch drops).  Random
mismatch is ngspice/idac_unit_size.py's job.
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
VDD = 1.8
I_BIAS = 50e-6          # 128 reference units x 390.6 nA
TRIM_NOMINAL = 64       # N = 64 + 64 = 128
C_SAR = 0.5e-12         # per side (spec: single-ended ramp into 0.5 pF)
#: extracted netlist to simulate instead of the schematic (--pex)
PEX = None
#: R/C process corner, e.g. 'res_high__cap_low' (--rc).  sky130's MOS corners
#: all load res_typical__cap_typical, so R and C only move if asked to.
RC = None
NG = "/foss/pdks/sky130A/libs.tech/ngspice"


def rc_lines() -> list[str]:
    if not RC:
        return []
    return [f'.include "{NG}/r+c/{RC}.spice"',
            f'.include "{NG}/r+c/{RC}__lin.spice"']
CBIAS = 1.30            # cascode gate: unit and cascode both saturated, ss..ff


def run(netlist: str) -> str:
    with tempfile.NamedTemporaryFile("w", suffix=".spice", delete=False) as fh:
        fh.write(netlist)
    r = subprocess.run(["ngspice", "-b", fh.name], capture_output=True,
                       text=True, timeout=3600)
    return r.stdout + r.stderr


def vals(text: str, name: str) -> list[float]:
    return [float(v) for v in re.findall(rf"^{name}\s*=\s*(\S+)", text, re.M)]


def header(corner: str, temp: float, trim: int) -> str:
    segs = dac.segments()
    conn = " ".join(dac.ports()).replace("ibias", "ib")
    lines = [f"* idac testbench, {corner} {temp} C",
             f".lib {LIB} {corner}",
             *rc_lines(),
             f".include {PEX}" if PEX else dac.netlist(),
             f".temp {temp}",
             # unselected trim units leave their drains floating; the shunt
             # gives them a defined voltage.  1 TOhm leaks ~1 pA against
             # 195 nA units.
             ".option rshunt=1e12",
             # only the outputs: the extracted netlist has ~2500 nodes, and
             # saving them all over a 20 us transient took 2.2 GB
             ".save v(outp) v(outn)",
             f"Vdd vdd 0 {VDD}",
             "Vss vss 0 0",
             f"Ib vdd ib {I_BIAS}",
             f"Xdut {conn} idac"]
    if dac.CASCODE:
        # ideal here; on chip it comes from the bias block, as the ramp
        # sink's does.  Sets the unit sources' Vds to ~0.5 V.
        lines.append(f"Vcb cbias 0 {CBIAS}")
    for k in range(len(dac.REF_TRIM)):
        lines.append(f"Vtrim{k} trim{k} 0 {VDD if trim >> k & 1 else 0}")
    for name, _ in segs:
        lines.append(f"Vd_{name} d_{name} 0 0")
        lines.append(f"Vdb_{name} db_{name} 0 {VDD}")
    return "\n".join(lines)


def code_bits(code: int) -> dict:
    """Switch state of every segment for one code: 1 = current to outp."""
    msb, lsb = code >> 6, code & 63
    on = {}
    for name, _ in dac.segments():
        if name.startswith("t"):
            on[name] = int(msb > int(name[1:]))
        else:
            on[name] = (lsb >> int(name[1:])) & 1
    return on


def pwl(levels: list[int], slot: float, mbb: float, edge: float = 0.1e-9) -> str:
    """PWL drive for one switch gate, one level per code slot.

    ``mbb`` delays every FALLING edge, so the switch turning on conducts
    before its partner turns off (make-before-break).  With both edges at
    once, mid-edge both switches are nearly off, the segment's source node
    collapses, and the kick couples through Cgd into the 12 pF ibias node.
    Times need ~9 significant digits: at 1024 x 400 ns a 1 ns offset is a few
    parts in 10^6, and printed to 5 digits point pairs collided.
    """
    pts, prev = [], None
    for code, b in enumerate(levels):
        if b != prev:
            t = code * slot + (mbb if prev == 1 and b == 0 else 0.0)
            if prev is not None:
                pts.append(f"{t:.9e} {prev * VDD}")
            pts.append(f"{t + edge:.9e} {b * VDD}")
            prev = b
    return " ".join(pts)


def boundary_codes() -> list[int]:
    """Endpoints plus two codes either side of every 64-code boundary: where
    the most switches move at once, and where the dynamic error concentrates.
    75 codes instead of 1024 -- what the extracted netlist can afford."""
    codes = {0, 1, 1022, 1023}
    for k in range(1, 16):
        codes.update(range(64 * k - 2, 64 * k + 3))
    return sorted(codes)


def static(corner: str, temp: float, trim: int, slot: float = 40e-9,
           mbb: float = 0.0, codes: list[int] | None = None) -> dict:
    """Every code, as one transient.

    Each code holds for ``slot`` and is sampled at the end of it.  A separate
    ``op`` per code restarts Newton from nothing and falls back to gmin
    stepping every time -- about a second each; stepping through the codes in
    a transient starts each one from the last solution instead.  ``slot`` of
    a few hundred ns gives the settled (static) transfer; 20 ns is the real
    update rate, where what is left unsettled is the dynamic error.
    """
    codes = list(range(1024)) if codes is None else codes
    body = header(corner, temp, trim)
    for name, _ in dac.segments():
        levels = [code_bits(c)[name] for c in codes]
        body = body.replace(f"Vd_{name} d_{name} 0 0",
                            f"Vd_{name} d_{name} 0 pwl({pwl(levels, slot, mbb)})")
        body = body.replace(
            f"Vdb_{name} db_{name} 0 {VDD}",
            f"Vdb_{name} db_{name} 0 pwl({pwl([1 - b for b in levels], slot, mbb)})")
    t_end = len(codes) * slot
    # The maximum step has to resolve the 0.1 ns switching edges.  At 2 ns the
    # integration error of every switching event survived to the sample point
    # and read as 0.35 LSB of "dynamic" DNL at 20 ns slots -- a circuit that a
    # 0.1 ns-step transient shows settled to 0.01 LSB within 1 ns.
    tmax = min(0.2e-9, slot / 100)
    text = run(body + f"""
.control
set wr_singlescale
set wr_vecnames
tran {tmax} {t_end} 0 {tmax}
wrdata /tmp/idac_static.txt v(outp) v(outn)
.endc
.end
""")
    try:
        lines = open("/tmp/idac_static.txt").read().split("\n")
    except OSError:
        raise RuntimeError(text[-3000:])
    data = [list(map(float, ln.split())) for ln in lines[1:] if ln.strip()]
    if not data:
        raise RuntimeError(text[-3000:])
    # sample each code 90 % into its slot
    vp, vn, j = [], [], 0
    for k in range(len(codes)):
        t = (k + 0.9) * slot
        while j + 1 < len(data) and data[j + 1][0] <= t:
            j += 1
        vp.append(data[j][1])
        vn.append(data[j][2])
    vd = [n - p for p, n in zip(vp, vn)]
    at = dict(zip(codes, vd))
    lsb = (at[1023] - at[0]) / 1023
    inl = [(at[c] - (at[0] + lsb * c)) / lsb for c in codes]
    dnl = [(at[c + 1] - at[c]) / lsb - 1 for c in codes if c + 1 in at]
    return {"lsb_v": lsb, "fs": at[1023] - at[0], "off": at[0],
            "codes": len(codes),
            "inl_max": max(map(abs, inl)), "dnl_max": max(map(abs, dnl), default=0.0),
            "cm_lo": min((p + n) / 2 for p, n in zip(vp, vn)),
            "cm_hi": max((p + n) / 2 for p, n in zip(vp, vn)),
            "v_lowest": min(min(vp), min(vn))}


def settling(corner: str, temp: float, trim: int) -> dict:
    """Full-scale step: all segments swing outp-side to outn-side at 2 ns."""
    segs = dac.segments()
    lines = [header(corner, temp, trim),
             f"Cp outp 0 {C_SAR}", f"Cn outn 0 {C_SAR}"]
    body = "\n".join(lines)
    for name, _ in segs:
        body = body.replace(f"Vd_{name} d_{name} 0 0",
                            f"Vd_{name} d_{name} 0 pulse(0 {VDD} 2n 20p 20p 1u 2u)")
        body = body.replace(f"Vdb_{name} db_{name} 0 {VDD}",
                            f"Vdb_{name} db_{name} 0 pulse({VDD} 0 2n 20p 20p 1u 2u)")
    text = run(body + """
.control
tran 5p 30n
let vd = v(outn) - v(outp)
let vf = vd[length(vd) - 1]
let v0 = vd[0]
let err = abs(vd - vf) / abs(vf - v0)
meas tran t01 when err=1e-3 fall=last
meas tran t63 when err=0.3679 fall=1
let settle_ns = (t01 - 2e-9) * 1e9
let tau_ns = (t63 - 2e-9) * 1e9
print settle_ns tau_ns vf v0
.endc
.end
""")
    out = {k: vals(text, k) for k in ("settle_ns", "tau_ns")}
    if not all(out.values()):
        raise RuntimeError(text[-2000:])
    return {k: v[0] for k, v in out.items()}


def headroom(corner: str, temp: float, trim: int) -> dict:
    """Operating point at mid-scale: the unit sources' Vds (cascode source) for
    a thermometer and the 1-unit segment, and the switch nodes above them."""
    nodes = ["xdut.c_t0", "xdut.c_b0", "xdut.s_t0", "xdut.s_b0"]
    body = header(corner, temp, trim).replace(".save v(outp) v(outn)", ".save all")
    text = run(body + f"""
.control
op
print {' '.join(f'v({n})' for n in nodes)}
.endc
.end
""")
    out = {}
    for n in nodes:
        v = vals(text, re.escape(f"v({n})"))
        out[n.split(".")[1]] = v[0] if v else float("nan")
    return out


def r_d(corner: str, temp: float) -> float:
    text = run(f"""* R_D
.lib {LIB} {corner}
{chr(10).join(rc_lines())}
.temp {temp}
XR a 0 0 {dac.RES} L={dac.RD_L}
Va a 0 1
.control
op
let r = 1 / (-i(Va))
print r
.endc
.end
""")
    return vals(text, "r")[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corner", default="tt")
    ap.add_argument("--temp", type=float, default=27.0)
    ap.add_argument("--trim", type=int, default=TRIM_NOMINAL)
    ap.add_argument("--skip-static", action="store_true")
    ap.add_argument("--pex", help="extracted netlist (layout/build/idac/idac_pex.spice)")
    ap.add_argument("--boundary", action="store_true",
                    help="75 codes around the 64-code boundaries instead of all 1024")
    args = ap.parse_args()
    global PEX
    PEX = args.pex

    res = {"corner": args.corner, "temp": args.temp, "trim": args.trim,
           "r_d_ohm": r_d(args.corner, args.temp)}
    print(f"{args.corner} {args.temp:g} C, trim {args.trim}")
    print(f"  R_D            {res['r_d_ohm']:.0f} Ohm")
    if not args.skip_static:
        # 20 ns per code: the real update rate, so what is left unsettled at
        # the sample point counts
        s = static(args.corner, args.temp, args.trim, slot=20e-9,
                   codes=boundary_codes() if args.boundary else None)
        res.update(s)
        print(f"  LSB            {s['lsb_v']*1e6:.1f} uV diff "
              f"(ideal 781.25 uV = 2 x I_u x R_D)")
        print(f"  full scale     {s['fs']*1e3:.1f} mV, offset {s['off']*1e3:.1f} mV")
        print(f"  INL / DNL      {s['inl_max']:.3f} / {s['dnl_max']:.3f} LSB "
              f"({s['codes']} codes at 20 ns each)")
        print(f"  common mode    {s['cm_lo']:.4f} .. {s['cm_hi']:.4f} V")
        print(f"  lowest output  {s['v_lowest']:.3f} V (ramp window floor 1.4 V)")
    t = settling(args.corner, args.temp, args.trim)
    res.update(t)
    print(f"  settling       {t['settle_ns']:.2f} ns to 0.1 %, "
          f"tau {t['tau_ns']:.2f} ns (budget dac_settle_tau <= 4 ns)")
    tag = "_pex" if args.pex else ""
    out = HERE / f"idac_tb_{args.corner}_{args.temp:g}{tag}.json"
    out.write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
