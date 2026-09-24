"""Fig. 5 -- best achievable in-band FDVPD noise floor versus its power.

At every power budget the split between ramp current and comparator power is
re-optimised, with the slew rate held at its nominal value (so ``C_SAR`` tracks
``I_R``).  The family of curves is a sweep over phase-detector resolution,
which is an ordered quantity -- hence the single-hue ramp rather than
categorical hues.

The marker is the prototype's own operating point.

    python scripts/fig05_noise_vs_power.py
"""

from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from fdvadpll import default_design, noise
from _style import INK_2, MUTED, SERIES, annotate, ordinal, save, title, use_style

BITS = (8, 10, 12, 14)
FRAC_BIT = 9


def main() -> None:
    use_style()
    design = default_design()
    powers = np.logspace(np.log10(0.2e-3), np.log10(20e-3), 120)
    colours = ordinal(len(BITS))

    fig, ax = plt.subplots(figsize=(7.6, 5.0))

    for bits, colour in zip(BITS, colours):
        L, info = noise.noise_vs_power(design, bits, powers,
                                       frac_bit=FRAC_BIT)
        ok = np.isfinite(L)
        ax.semilogx(powers[ok] * 1e3, L[ok], color=colour,
                    label=f"{bits} bit PD")
        # direct label at the right-hand end -- relief for the low-contrast steps
        ax.annotate(f"{bits} b", xy=(powers[ok][-1] * 1e3, L[ok][-1]),
                    xytext=(6, 0), textcoords="offset points", fontsize=8.5,
                    color=INK_2, va="center")

    # the prototype: FDVPD power from the paper's Fig. 17 share of 9.2 mW
    p = design.power
    p_fdvpd = (p.share_dac_ramp + p.share_adc) * p.p_total_measured
    floor = 10 * np.log10(noise.pd_noise_floor(design, FRAC_BIT).pd_total)
    # upper right of the marker is the only clear space between the curves
    annotate(ax, p_fdvpd * 1e3, floor,
             f"  prototype\n  {p_fdvpd*1e3:.2f} mW, {floor:.0f} dBc/Hz",
             SERIES[1], dx=8, dy=8, ha="left", va="bottom")

    title(ax, "In-band FDVPD noise floor versus phase-detector power",
          f"slew rate fixed at {design.fdvpd.slew_rate*1e-9:.1f} mV/ps; "
          f"I_R / comparator split re-optimised at every budget")
    ax.set_xlabel("FDVPD power [mW]")
    ax.set_ylabel("in-band L(f) [dBc/Hz]")
    ax.set_xlim(0.2, 26)
    ax.set_xticks([0.2, 0.5, 1, 2, 5, 10, 20])
    ax.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10", "20"])
    ax.minorticks_off()
    # direct labels already sit at the right-hand ends, so the legend goes
    # where nothing else is
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "fig05_noise_vs_power.png")

    print("\n  optimum device sizing at a few budgets:")
    probe = np.array([0.5e-3, 1e-3, 2e-3, 5e-3])
    L, info = noise.noise_vs_power(design, 12, probe, frac_bit=FRAC_BIT)
    for i, budget in enumerate(probe):
        print(f"    {budget*1e3:4.1f} mW -> {L[i]:7.2f} dBc/Hz   "
              f"I_R {info['i_ramp'][i]*1e6:6.1f} uA   "
              f"C_SAR {info['c_sar'][i]*1e15:6.0f} fF   "
              f"v_TN {info['v_tn'][i]*1e6:5.1f} uV")


if __name__ == "__main__":
    main()
