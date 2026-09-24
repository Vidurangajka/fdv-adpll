"""The background calibration loops (the ADPLL extension).

Three panels:

1. Acquisition -- the counter-based FLL pulls a 30 MHz initial error in, then
   hands over to the phase detector once the ADC stops saturating.
2. PD gain calibration -- a 1 % DAC-to-ramp mismatch turns the fractional
   accumulator's sawtooth into output phase modulation; the LMS trims it out.
3. The spur before and after, measured rather than asserted.

    python scripts/calibration_demo.py
"""

from __future__ import annotations

import sys
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from fdvadpll import FdvPll, fractional_design, measured_fit_design
from _style import AXIS, INK_2, MUTED, SERIES, save, title, use_style

N_CYCLES = 1 << 16
FRAC_BIT = 5
GAIN_ERR = 0.01


def main() -> None:
    use_style()
    base = measured_fit_design()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    # ---- 1: acquisition --------------------------------------------------
    print("  acquisition ...")
    acq = FdvPll(base, seed=2, f_init_offset=-30e6).run(1 << 15, discard=0)
    t_us = np.arange(len(acq.f_dco)) / base.f_ref * 1e6
    ax = axes[0]
    ax.plot(t_us, (acq.f_dco - base.f_ckv) / 1e6, color=SERIES[0],
            linewidth=1.2)
    ax.axhline(0.0, color=AXIS, linewidth=0.9)
    if acq.fll_disabled_at is not None:
        t_off = acq.fll_disabled_at / base.f_ref * 1e6
        ax.axvline(t_off, color=MUTED, linewidth=0.9, linestyle=(0, (2, 3)))
        ax.annotate(f"FLL off\n{t_off:.1f} us", xy=(t_off, 0), xytext=(6, -28),
                    textcoords="offset points", fontsize=8.5, color=INK_2)
    title(ax, "Frequency acquisition", "30 MHz initial error")
    ax.set_xlabel("time [us]")
    ax.set_ylabel("f_DCO - f_CKV [MHz]")
    ax.set_xlim(0, 8)          # the whole transient is over inside 5 us

    # ---- 2: gain calibration --------------------------------------------
    print("  gain calibration ...")
    d = fractional_design(FRAC_BIT, base)
    d_bad = replace(d, fdvpd=replace(d.fdvpd, dac_gain_err=GAIN_ERR))
    cal = FdvPll(d_bad, seed=5, calibrate=("gain",),
                 cal_mu={"gain": 5e-3}).run(N_CYCLES)
    ax = axes[1]
    t_us = np.arange(len(cal.gain_history)) / base.f_ref * 1e6
    ax.plot(t_us, cal.gain_history, color=SERIES[1], linewidth=1.4,
            label="estimate")
    target = 1.0 / (1.0 + GAIN_ERR)
    ax.axhline(target, color=AXIS, linewidth=1.0, linestyle=(0, (4, 3)))
    ax.annotate(f"ideal {target:.4f}", xy=(t_us[-1], target), xytext=(-4, 6),
                textcoords="offset points", ha="right", fontsize=8.5,
                color=INK_2)
    title(ax, "PD gain calibration", f"{GAIN_ERR*100:.0f} % DAC-to-ramp mismatch")
    ax.set_xlabel("time [us]")
    ax.set_ylabel("gain correction  g")

    # ---- 3: the spur, before and after ----------------------------------
    print("  spur comparison ...")
    ideal = FdvPll(d, seed=5).run(N_CYCLES)
    off = FdvPll(d_bad, seed=5).run(N_CYCLES)
    ax = axes[2]
    f_frac = abs(d.fcw - round(d.fcw)) * base.f_ref
    rows = [("ideal I-DAC", ideal, SERIES[0]),
            ("1 % mismatch", off, SERIES[1]),
            ("mismatch + cal.", cal, SERIES[2])]
    for i, (label, res, colour) in enumerate(rows):
        f, dbc = res.spurs()
        m = (f > 1e5) & (f < 2e7)
        ax.semilogx(f[m], dbc[m], color=colour, linewidth=0.9, alpha=0.8,
                    label=f"{label}  ({res.fractional_spur()[1]:.0f} dBc)")
    ax.axvline(f_frac, color=MUTED, linewidth=0.9, linestyle=(0, (2, 3)))
    ax.annotate(f"f_frac\n{f_frac/1e6:.1f} MHz", xy=(f_frac, -40),
                xytext=(6, 0), textcoords="offset points", fontsize=8.5,
                color=INK_2)
    title(ax, "Spectrum at the fractional frequency",
          f"FCW_pd = {d.fcw_pd:.5f}")
    ax.set_xlabel("offset frequency [Hz]")
    ax.set_ylabel("level [dBc]")
    ax.set_ylim(-140, -30)
    ax.legend(loc="lower left")

    fig.tight_layout()
    save(fig, "calibration_demo.png")

    print(f"\n    ideal I-DAC     {ideal.fractional_spur()[1]:7.1f} dBc")
    print(f"    1 % mismatch    {off.fractional_spur()[1]:7.1f} dBc")
    print(f"    + calibration   {cal.fractional_spur()[1]:7.1f} dBc")
    print(f"    gain estimate   {cal.gain_history[-1]:.5f} "
          f"(ideal {target:.5f})")


if __name__ == "__main__":
    main()
