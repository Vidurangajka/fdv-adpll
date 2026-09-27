#!/usr/bin/env python3
"""Current DAC across PVT and R/C corners (schematic).

    .\\osic.ps1 python3 ngspice/idac_corners.py        # ~40 min

MOS corners ss / tt / ff at -40 / 27 / 85 C with typical R and C, then the four
R/C corners at tt, 27 C -- sky130's MOS corners never move R or C, and R_D * C
is what the trim exists for.  Each corner: R_D, LSB, INL / DNL on the 79
boundary codes at 20 ns per code, settling into 0.5 pF, and the unit sources'
headroom.

``trim needed`` is the code that puts the DAC-to-ramp gain back on target:
V_DM has to span SR * T_PD, i.e. LSB = I_sink T_PD / (1024 C_SAR), so with the
corner's C_SAR (ngspice/rc_spread.sh) the unit current has to scale by
LSB_target / LSB_measured, and I_u = I_bias / (64 + trim).  It has to land in
0 .. 127 at every corner, with room for temperature on top.

Writes ngspice/idac_corners.json.
"""

from __future__ import annotations

import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import idac_tb as tb  # noqa: E402

#: C_SAR (0.5 pF drawn) per R/C corner at 27 C: rc_spread.sh's 1 pF MIM, halved
C_SAR = {rc: 0.5 * c for rc, c in {
    None: 0.9825e-12, "res_high__cap_high": 1.1036e-12,
    "res_high__cap_low": 0.8562e-12, "res_low__cap_high": 1.1146e-12,
    "res_low__cap_low": 0.8648e-12}.items()}
C_NOM = 0.5e-12
LSB_NOM = 781.25e-6             # 200 uA * 2 ns / (1024 * 0.5 pF)
TRIM = 64                       # N = 128, the nominal


def one(corner, temp, rc):
    tb.RC = rc
    r = {"corner": corner, "temp": temp, "rc": rc or "typical",
         "r_d_ohm": tb.r_d(corner, temp)}
    r.update(tb.static(corner, temp, TRIM, slot=20e-9, codes=tb.boundary_codes()))
    r.update(tb.settling(corner, temp, TRIM))
    r.update(tb.headroom(corner, temp, TRIM))
    target = LSB_NOM * C_NOM / C_SAR[rc]
    n_needed = (64 + TRIM) * r["lsb_v"] / target
    r["trim_needed"] = n_needed - 64
    return r


def main():
    runs = [(c, t, None) for c in ("ss", "tt", "ff") for t in (-40, 27, 85)]
    runs += [("tt", 27, rc) for rc in C_SAR if rc]
    rows = []
    print(f"{'corner':24s} {'R_D':>6} {'LSB uV':>7} {'INL':>6} {'DNL':>6} "
          f"{'tau ns':>6} {'Vout min':>8} {'Vds t0':>7} {'Vds b0':>7} {'trim':>6}")
    for c, t, rc in runs:
        r = one(c, t, rc)
        rows.append(r)
        name = f"{c} {t:+d}C" + (f" {rc}" if rc else "")
        print(f"{name:24s} {r['r_d_ohm']:6.0f} {r['lsb_v']*1e6:7.1f} "
              f"{r['inl_max']:6.3f} {r['dnl_max']:6.3f} {r['tau_ns']:6.2f} "
              f"{r['v_lowest']:8.3f} {r['c_t0']:7.3f} {r['c_b0']:7.3f} "
              f"{r['trim_needed']:6.1f}", flush=True)
    (HERE / "idac_corners.json").write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
