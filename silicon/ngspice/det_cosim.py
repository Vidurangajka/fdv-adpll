#!/usr/bin/env python3
"""The whole detector, mixed-signal: the control RTL runs the analog blocks.

    .\\osic.ps1 python3 ngspice/det_cosim.py [--cycles 24] [--drift 0.6]

ngspice's d_cosim bridge runs digital/fdvpd_ctrl.v (Icarus, -DSIM for the
clock gates and delay lines) against transistor netlists of the DAC
(idac_netlist), the ramp sink, and the SAR front end (fe_netlist: S1, the
ramp's steering switch, the CDAC and the comparator).  Nothing between them is idealised but the bias:
ideal ibias / cbias / ramp / ramp-switch currents, and an ideal reference
ladder.

CKV runs at FCW = 2 x (9 + frac), frac = 57819 / 65536 by default: the DAC
code steps ~121 a cycle and wraps every ~8.5, so it is not collinear with the
phase, which starts `phase` ps past the ideal model's lock point and drifts
by `drift` ps per reference cycle -- the residue sweeps the SAR's range while
the code walks through the DAC.  Each cycle the log gives the DAC
code, the S2 window dt and the SAR code; the testbench also reads vn - vp at
the first comparison.  Two checks:

  SAR       the code against an ideal 7 b quantisation of that residue
            (LSB 195.3 uV, offset binary) -- the CDAC, ladder and comparator
  detector  code = a * dac + b * dt + c over the unrailed cycles: 1/b is the
            detector's time LSB, -a/b / (T_PD/1024) the DAC-to-ramp gain
            match (1 at trim), and the fit residual its linearity and noise
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
SILICON = HERE.parent
LAYOUT = SILICON / "layout"
DIGITAL = SILICON / "digital"
sys.path.insert(0, str(LAYOUT))
import fe_netlist as fe  # noqa: E402
import idac_netlist as dac  # noqa: E402

LIB = "/foss/pdks/sky130A/libs.tech/ngspice/sky130.lib.spice"
VDD = 1.8
T_REF = 20e-9
T_PD = 2e-9
FCW_INT = 9
#: FCW / 2 = 9 + frac: 1 - frac = 0.1178, so the DAC code steps ~121 a
#: cycle and wraps every ~8.5 -- uncorrelated with the linear phase drift
FCW_FRAC = 57819
#: the DAC-to-ramp gain match at tt: 64, nominal (70 compensated the old PMOS
#: ramp switch, which slowed the ramp ~6 %)
TRIM = 64
LSB_ADC = 195.3125e-6
VRM = 0.75
#: ladder step: vrp20 - vrm is eight of them, one LSB-weight unit is two
V_STEP = LSB_ADC * fe.C_TOT / (2 * fe.C_UNIT) / 2
TAPS = {"vrp20": 8, "vrp10": 4, "vrp5": 2, "vrm": 0, "vrn5": -2, "vrn10": -4,
        "vrn20": -8}

IN_PORTS = (["rst_n", "ref_clk", "ckv"] + [f"fi{k}" for k in range(7, -1, -1)]
            + [f"ff{k}" for k in range(15, -1, -1)]
            + [f"ti{k}" for k in range(6, -1, -1)] + ["cmp_p", "cmp_n"])
DSW = [f"b{k}" for k in range(6)] + [f"t{k}" for k in range(14, -1, -1)]
OUT_PORTS = (["ckvd"] + [f"d_{s}" for s in DSW] + [f"db_{s}" for s in DSW]
             + [f"trim{k}" for k in range(6, -1, -1)]
             + ["s1_n", "s2_n", "cmp_clk"]
             + [f"bpos{k}" for k in range(5, -1, -1)]
             + [f"bneg{k}" for k in range(5, -1, -1)]
             + [f"bmid{k}" for k in range(5, -1, -1)]
             + [f"code{k}" for k in range(6, -1, -1)]
             + ["code_valid", "railed"] + [f"cref{k}" for k in range(7, -1, -1)])
#: outputs bridged to the analog side; dswb and trim come another way (below)
ANALOG_OUT = [p for p in OUT_PORTS if not p.startswith(("db_", "trim"))]
#: what --debug writes to dbg.txt: the analog nodes and every DAC input
DEBUG = (["ckv", "ref_clk", "outp", "outn", "rs", "ib", "cbias", "s1_n", "s2_n", "cmp_clk", "vg",
          "cmp_p", "cmp_n", "ckvd"] + [f"bpos{k}" for k in range(6)]
         + [f"bneg{k}" for k in range(6)] + [f"d_{s}" for s in DSW]
         + [f"db_{s}" for s in DSW] + [f"trim{k}" for k in range(7)])


def clocks(drift_ps: float, phase_ps: float, frac: int) -> tuple[float, float, float]:
    """(t_ckv0, t_ckv, t_ref0), s: the first CKV rise, its period, the first REF."""
    fcw = FCW_INT + frac / 65536
    t_ckv = (T_REF + drift_ps * 1e-12) / fcw / 2
    t_ckv0 = 1e-9
    # the gate passes ckv edge 19 first; REF 300 ps (code 0) before it, less
    # 18 ps, and phase_ps earlier still
    t_ref0 = t_ckv0 + 18 * t_ckv - 300e-12 + 18e-12 - phase_ps * 1e-12
    return t_ckv0, t_ckv, t_ref0


def clock_defines(drift_ps: float, phase_ps: float, frac: int) -> list[str]:
    """iverilog -D flags that make cosim_top generate CKV and REF itself, at
    the times the sources below cross 0.9 V (10 ps into a 20 ps edge, plus the
    bridge's 1 ps) -- see cosim_top.v for why they are not bridged."""
    t_ckv0, t_ckv, t_ref0 = clocks(drift_ps, phase_ps, frac)
    ps = lambda t: f"{t * 1e12:.6f}"
    return [f"-DTCKV0={ps(t_ckv0 + 11e-12)}", f"-DTCKV={ps(t_ckv)}",
            f"-DTREF0={ps(t_ref0 + 11e-12)}", f"-DTREF={ps(T_REF)}"]


def netlist(cycles: int, drift_ps: float, workdir: pathlib.Path, phase_ps: float = 0.0,
            debug: bool = False, frac: int = FCW_FRAC,
            trim: int = TRIM) -> tuple[str, float]:
    t_ckv0, t_ckv, t_ref0 = clocks(drift_ps, phase_ps, frac)
    t_end = t_ref0 + (cycles + 1) * T_REF
    L = [f"* detector co-simulation, {cycles} cycles", f".lib {LIB} tt",
         dac.netlist(), (LAYOUT / "ramp_sink_ref.spice").read_text(),
         fe.netlist(), ".option rshunt=1e12",
         f"Vdd vdd 0 {VDD}", "Vss vss 0 0",
         "Ib vdd ib 50u", "Vcb cbias 0 1.30",
         "Irr vdd nbias 200u", "Irc vdd rcbias 200u",
         # the ramp switch's gate bias: 20 uA into the front end's replica
         "Isw vdd vg 20u"]
    for tap, n in TAPS.items():
        L.append(f"V{tap} {tap} 0 {VRM + n * V_STEP:.6f}")
    # clocks and static inputs, as voltages into adc bridges
    L += [f"Vckv ckv 0 pulse(0 {VDD} {t_ckv0:.6e} 20p 20p {t_ckv/2 - 20e-12:.6e} {t_ckv:.9e})",
          f"Vref ref_clk 0 pulse(0 {VDD} {t_ref0:.6e} 20p 20p {T_REF/2 - 20e-12:.6e} {T_REF:.6e})",
          f"Vrst rst_n 0 pwl(0 {VDD} 10p {VDD} 30p 0 480p 0 500p {VDD})"]
    for k in range(8):
        L.append(f"Vfi{k} fi{k} 0 {VDD if FCW_INT >> k & 1 else 0}")
    for k in range(16):
        L.append(f"Vff{k} ff{k} 0 {VDD if frac >> k & 1 else 0}")
    for k in range(7):
        L.append(f"Vti{k} ti{k} 0 {VDD if trim >> k & 1 else 0}")
        # trim passes straight through the block: drive the DAC from DC too
        L.append(f"Vtrim{k} trim{k} 0 {VDD if trim >> k & 1 else 0}")
    dins = [f"x_{p}" for p in IN_PORTS]
    douts = [f"y_{p}" for p in OUT_PORTS]
    L += [f"aadc [{' '.join(IN_PORTS)}] [{' '.join(dins)}] a2d",
          ".model a2d adc_bridge in_low=0.9 in_high=0.9 rise_delay=1e-12 fall_delay=1e-12",
          f"actl [{' '.join(dins)}] [{' '.join(douts)}] null ctl",
          f'.model ctl d_cosim simulation="ivlng" sim_args=["{workdir}/ctl.vvp"] delay=1e-12',
          f"adac [{' '.join(f'y_{p}' for p in ANALOG_OUT)}] [{' '.join(ANALOG_OUT)}] d2a",
          f".model d2a dac_bridge out_low=0 out_high={VDD} t_rise=30p t_fall=30p",
          # the dump gates from dsw, inverted: dswb is ~dsw in the RTL, and a
          # bridge's output is 0 V at the operating point -- through dswb's own
          # bridge every unit would start cut off, and the DAC's bias node,
          # ~1300 unit gates on a 50 uA diode, takes microseconds to recover
          f"adacb [{' '.join(f'y_d_{s}' for s in DSW)}] [{' '.join(f'db_{s}' for s in DSW)}] d2ab",
          f".model d2ab dac_bridge out_low={VDD} out_high=0 t_rise=30p t_fall=30p"]
    # the analog blocks
    L.append(f"Xdac {' '.join(dac.ports()).replace('ibias', 'ib')} idac")
    L.append("Xramp vss rs nbias rcbias ramp_sink")
    L.append(f"Xfe {' '.join(fe.ports())} sarfe")
    L += [".control", "set wr_singlescale", "set wr_vecnames",
          f"tran 10p {t_end:.6e} 0 20p",
          f"wrdata {workdir}/vpn.txt v(vp) v(vn)",
          *([f"wrdata {workdir}/dbg.txt " + " ".join(f"v({n})" for n in DEBUG)]
            if debug else []),
          ".endc", ".end"]
    return "\n".join(L) + "\n", t_end


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cycles", type=int, default=24)
    ap.add_argument("--drift", type=float, default=1.9, help="ps per REF cycle")
    ap.add_argument("--frac", type=int, default=FCW_FRAC, help="FCW/2 fraction, 2^-16")
    ap.add_argument("--trim", type=int, default=TRIM)
    ap.add_argument("--phase", type=float, default=70.0,
                    help="ps the first window is longer than the ideal model's lock point")
    ap.add_argument("--debug", action="store_true",
                    help="also write every analog node and DAC input to dbg.txt")
    ap.add_argument("--work", type=pathlib.Path, default=SILICON / "build" / "det_cosim")
    args = ap.parse_args()
    w = args.work
    w.mkdir(parents=True, exist_ok=True)
    subprocess.run(["iverilog", "-g2012", "-DSIM",
                    *clock_defines(args.drift, args.phase, args.frac),
                    "-o", str(w / "ctl.vvp"),
                    str(DIGITAL / "cosim_top.v"), str(DIGITAL / "fdvpd_ctrl.v")],
                   check=True)
    text, _ = netlist(args.cycles, args.drift, w, args.phase, args.debug, args.frac,
                      args.trim)
    (w / "det.cir").write_text(text)
    log = w / "cosim_log.txt"
    log.unlink(missing_ok=True)
    out = subprocess.run(["ngspice", "-b", "det.cir"], cwd=w, capture_output=True,
                         text=True, timeout=6 * 3600)
    if not log.exists() or not (w / "vpn.txt").exists():
        raise RuntimeError(out.stdout[-3000:] + out.stderr[-3000:])

    rows = np.loadtxt(log, ndmin=2)
    wave = np.loadtxt(w / "vpn.txt", skiprows=1)
    t, vp, vn = wave[:, 0], wave[:, 1], wave[:, 2]
    n, dcode, dt, tc, code, railed = (rows[:, 0], rows[:, 1], rows[:, 3] * 1e-12,
                                      rows[:, 4] * 1e-12, rows[:, 5], rows[:, 6])
    vd = np.interp(tc - 20e-12, t, vn - vp)
    ideal = np.clip(np.floor(vd / LSB_ADC + 64), 0, 127)
    keep = n > 2                              # the first cycles start from reset
    sar_err = (code - ideal)[keep]
    print(f"cycles {int(n[-1])} logged, {int(keep.sum())} checked")
    print(f"{'n':>3} {'dac':>5} {'dt ps':>8} {'vd mV':>8} {'code':>5} {'ideal':>6}")
    for i in range(len(n)):
        print(f"{int(n[i]):3d} {int(dcode[i]):5d} {dt[i]*1e12:8.1f} "
              f"{vd[i]*1e3:8.3f} {int(code[i]):5d} {int(ideal[i]):6d}")
    print(f"SAR       code - ideal: {np.unique(sar_err, return_counts=True)}")
    ok = keep & (code > 0) & (code < 127)
    res = {"cycles": rows.tolist(), "vd_v": vd.tolist(), "ideal": ideal.tolist()}
    if ok.sum() >= 5:
        A = np.column_stack([dcode[ok], dt[ok] * 1e12, np.ones(ok.sum())])
        coef, *_ = np.linalg.lstsq(A, code[ok], rcond=None)
        err = code[ok] - A @ coef
        lsb_t = -1 / coef[1]
        match = -coef[0] / coef[1] / (T_PD / 1024 * 1e12)
        print(f"detector  {lsb_t*1e3:.0f} fs per code, DAC/ramp gain match {match:.4f}, "
              f"fit residual rms {np.std(err):.2f} LSB, max {np.max(np.abs(err)):.2f} "
              f"({ok.sum()} unrailed cycles)")
        res["fit"] = {"code_per_dac": coef[0], "code_per_ps": coef[1],
                      "offset": coef[2], "lsb_fs": lsb_t * 1e3, "gain_match": match,
                      "residual_rms_lsb": float(np.std(err))}
    (HERE / "det_cosim.json").write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
