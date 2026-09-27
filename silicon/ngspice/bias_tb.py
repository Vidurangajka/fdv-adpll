#!/usr/bin/env python3
"""Characterise the bias block (layout/bias_ref.spice) with its real loads.

    .\\osic.ps1 python3 ngspice/bias_tb.py            # corners + start-up
    .\\osic.ps1 python3 ngspice/bias_tb.py --tune     # size the free devices

Loads: the ramp sink's own netlist on iramp / icbr (outp held at 1.3 V, the
middle of its swing), the DAC's 128-unit reference diode on idac, nothing on
dcbias (it drives gates), and the ladder as drawn.  Reports each output
current against its target, dcbias against the 1.30 V the DAC's headroom
sweep settled on, the ladder span against 16 CDAC steps (42.0 mV), and whether the core
starts from a slow VDD ramp.

The absolute current sets the ramp slope, SR = I_R / C_SAR, and so the
detector's gain; the spread across corners is what the loop has to live with
(the DAC-to-ramp ratio does not move: both mirror from here).
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
LAYOUT = HERE.parent / "layout"
LIB = "/foss/pdks/sky130A/libs.tech/ngspice/sky130.lib.spice"
NG = "/foss/pdks/sky130A/libs.tech/ngspice"
sys.path.insert(0, str(LAYOUT))
import fe_netlist as fe  # noqa: E402

#: the free sizes, substituted as literal numbers -- sky130's device
#: subcircuits break on {expressions}
SIZES = {"CBM": "6", "CBL": "1", "RLAD": "4.81"}
TARGET = {"iramp": 200e-6, "icbr": 200e-6, "idac": 50e-6}
DCBIAS_TARGET = 1.30
#: 16 ladder units; one unit is half an LSB weight's step on a unit MIM
LADDER_SPAN = 16 * 195.3125e-6 * fe.C_TOT / (4 * fe.C_UNIT)


def netlist(sizes: dict) -> str:
    s = (LAYOUT / "bias_ref.spice").read_text()
    for k, v in sizes.items():
        s = re.sub(rf"\b{k}\b", v, s)
    return s


def run(text: str) -> str:
    with tempfile.NamedTemporaryFile("w", suffix=".spice", delete=False) as fh:
        fh.write(text)
    return subprocess.run(["ngspice", "-b", fh.name], capture_output=True,
                          text=True, timeout=1800).stdout


def vals(text, name):
    return [float(v) for v in re.findall(rf"^{re.escape(name)}\s*=\s*(\S+)", text, re.M)]


def bench(corner: str, temp: float, sizes: dict, rc: str | None = None) -> dict:
    rcl = (f'.include "{NG}/r+c/{rc}.spice"\n.include "{NG}/r+c/{rc}__lin.spice"'
           if rc else "")
    text = run(f"""* bias bench
.lib {LIB} {corner}
{rcl}
.temp {temp}
{netlist(sizes)}
.include {LAYOUT / "ramp_sink_ref.spice"}
Vdd vdd 0 1.8
Vss vss 0 0
Xb vdd vss ir icb id dcbias vrp20 vrp10 vrp5 vrm vrn5 vrn10 vrn20 isw bias
* isw into a 1.6 V stand-in for the front end's vg
Visw isw 0 1.6
* ammeters into the real loads
Va1 ir nb 0
Va2 icb cb 0
Va3 id nd 0
Xr vss vo nb cb ramp_sink
Vo vo 0 1.3
XDd nd nd 0 0 sky130_fd_pr__nfet_01v8 W=0.84 L=24 m=128
.control
op
print i(va1) i(va2) i(va3) v(dcbias) v(vrp20) v(vrn20) v(vrm) v(vrp5)
.endc
.end
""")
    keys = ["i(va1)", "i(va2)", "i(va3)", "v(dcbias)", "v(vrp20)", "v(vrn20)",
            "v(vrm)", "v(vrp5)"]
    out = {k: vals(text, k) for k in keys}
    if not all(out.values()):
        raise RuntimeError(text[-2500:])
    o = {k: v[0] for k, v in out.items()}
    return {"corner": corner, "temp": temp, "rc": rc or "typical",
            "iramp": o["i(va1)"], "icbr": o["i(va2)"], "idac": o["i(va3)"],
            "dcbias": o["v(dcbias)"], "ladder_span": o["v(vrp20)"] - o["v(vrn20)"],
            "ladder_mid": o["v(vrm)"],
            "ladder_lsb_tap": o["v(vrp5)"] - o["v(vrm)"]}


def startup(sizes: dict, corner: str = "tt", temp: float = 27) -> dict:
    """VDD from 0 to 1.8 V over 10 us, then hold: does the core wake up?"""
    text = run(f"""* start-up
.lib {LIB} {corner}
.temp {temp}
{netlist(sizes)}
.include {LAYOUT / "ramp_sink_ref.spice"}
Vdd vdd 0 pwl(0 0 10u 1.8 40u 1.8)
Vss vss 0 0
Xb vdd vss ir icb id dcbias vrp20 vrp10 vrp5 vrm vrn5 vrn10 vrn20 isw bias
Visw isw 0 1.6
Va1 ir nb 0
Xr vss vo nb cb ramp_sink
Rcb icb cb 1
Vo vo 0 1.3
XDd id id 0 0 sky130_fd_pr__nfet_01v8 W=0.84 L=24 m=128
.control
tran 20n 40u
meas tran i_end find i(va1) at=39u
meas tran t_up when i(va1)=100u rise=1
.endc
.end
""")
    i_end = vals(text, "i_end")
    t_up = vals(text, "t_up")
    return {"iramp_end": i_end[0] if i_end else float("nan"),
            "t_up_s": t_up[0] if t_up else float("nan")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tune", action="store_true")
    args = ap.parse_args()
    if args.tune:
        for rl in ("4.6", "4.8", "5.0"):
            for cbm, cbl in (("4", "1"), ("6", "1"), ("8", "1")):
                sz = {"CBM": cbm, "CBL": cbl, "RLAD": rl}
                r = bench("tt", 27, sz)
                print(f"RLAD {rl} CB {cbm} x 5/{cbl}: iramp {r['iramp']*1e6:6.1f} uA  "
                      f"dcbias {r['dcbias']:.3f}  ladder {r['ladder_span']*1e3:.2f} mV "
                      f"mid {r['ladder_mid']:.3f}", flush=True)
        return
    rows = []
    runs = [(c, t, None) for c in ("ss", "tt", "ff") for t in (-40, 27, 85)]
    runs += [("tt", 27, rc) for rc in ("res_high__cap_high", "res_high__cap_low",
                                       "res_low__cap_high", "res_low__cap_low")]
    print(f"{'corner':28s} {'iramp':>7} {'icbr':>7} {'idac':>6} {'dcbias':>7} "
          f"{'ladder':>7} {'mid':>6}")
    for c, t, rc in runs:
        r = bench(c, t, SIZES, rc)
        rows.append(r)
        name = f"{c} {t:+d}C" + (f" {rc}" if rc else "")
        print(f"{name:28s} {r['iramp']*1e6:7.1f} {r['icbr']*1e6:7.1f} "
              f"{r['idac']*1e6:6.1f} {r['dcbias']:7.3f} "
              f"{r['ladder_span']*1e3:7.2f} {r['ladder_mid']:6.3f}", flush=True)
    s = {}
    for c in ("ss", "tt", "ff"):
        for t in (-40, 27, 85):
            s[f"{c}{t:+d}"] = st = startup(SIZES, c, t)
            print(f"start-up {c} {t:+d}C: iramp {st['iramp_end']*1e6:.1f} uA at 39 us, "
                  f"reached 100 uA at {st['t_up_s']*1e6:.2f} us", flush=True)
    (HERE / "bias_tb.json").write_text(json.dumps({"sizes": SIZES, "corners": rows,
                                                    "startup": s}, indent=2))


if __name__ == "__main__":
    main()
