"""Fig. 17 and Table I -- where the 9.2 mW goes, and what it buys.

Left:  the measured power breakdown.  The oscillator and its divide-by-2 take
       nearly three quarters of the budget, which is the point of the paper --
       the phase detector that replaces a TDC costs almost nothing.
Right: the jitter-power trade, with the prototype's operating point marked.
       Contours of constant FoM show how much of the result is the detector and
       how much is simply the power spent.

    python scripts/fig17_power_and_fom.py
"""

from __future__ import annotations

import math
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from fdvadpll import measured_fit_design, noise
from _style import (AXIS, INK, INK_2, MUTED, SERIES, annotate, save, title,
                    use_style)

FOM_CONTOURS = (-240, -245, -250, -255)


def main() -> None:
    use_style()
    d = measured_fit_design()
    p = d.power

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # ---- left: the power breakdown --------------------------------------
    blocks = [
        ("DCO + div-2", p.share_dco_div2, SERIES[0]),
        ("I-DAC + ramp", p.share_dac_ramp, SERIES[1]),
        ("SAR ADC", p.share_adc, SERIES[2]),
        ("digital loop filter", p.share_dlf, SERIES[3]),
        ("everything else", p.share_rest, SERIES[4]),
    ]
    y = np.arange(len(blocks))
    mw = [s * p.p_total_measured * 1e3 for _, s, _ in blocks]
    ax.barh(y, mw, color=[c for _, _, c in blocks], height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels([n for n, _, _ in blocks])
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    for yi, (name, share, _), v in zip(y, blocks, mw):
        ax.annotate(f"{v:.2f} mW   {share*100:.0f} %", xy=(v, yi),
                    xytext=(6, 0), textcoords="offset points", va="center",
                    fontsize=8.5, color=INK_2)
    pd_share = p.share_dac_ramp + p.share_adc
    title(ax, "Measured power breakdown",
          f"{p.p_total_measured*1e3:.1f} mW total  |  the whole phase detector "
          f"is {pd_share*100:.0f} % of it "
          f"({pd_share*p.p_total_measured*1e3:.2f} mW)")
    ax.set_xlabel("power [mW]")
    ax.set_xlim(0, max(mw) * 1.45)

    # ---- right: jitter versus power, with FoM contours -------------------
    powers = np.logspace(math.log10(0.3e-3), math.log10(60e-3), 200)
    for fom_db in FOM_CONTOURS:
        # FoM = 20log10(sigma) + 10log10(P/1mW)  ->  sigma(P)
        sigma = 10 ** ((fom_db - 10 * np.log10(powers / 1e-3)) / 20.0)
        ax2.loglog(powers * 1e3, sigma * 1e15, color=AXIS, linewidth=0.9,
                   linestyle=(0, (3, 3)), zorder=1)
        # label at the left end: the right-hand ends run into the axis
        ax2.annotate(f"{fom_db} dB", xy=(powers[0] * 1e3, sigma[0] * 1e15),
                     xytext=(2, 4), textcoords="offset points", fontsize=7.5,
                     color=MUTED, ha="left")

    f = np.logspace(3, math.log10(40e6), 4000)
    j = noise.integrated_jitter(f, noise.output_phase_noise(d, f, None),
                                d.f_ckv)
    annotate(ax2, p.p_total_measured * 1e3, j * 1e15,
             f"  this design\n  {j*1e15:.0f} fs at "
             f"{p.p_total_measured*1e3:.1f} mW\n  FoM "
             f"{noise.fom(j, p.p_total_measured):.1f} dB",
             SERIES[0], dy=-6, va="top")

    title(ax2, "Jitter-power trade", "dashed lines are constant FoM")
    ax2.set_xlabel("total power [mW]")
    ax2.set_ylabel("rms jitter [fs]")
    ax2.set_xlim(0.3, 60)
    ax2.set_ylim(20, 2000)
    ax2.set_xticks([0.5, 1, 5, 10, 50])
    ax2.set_xticklabels(["0.5", "1", "5", "10", "50"])
    ax2.set_yticks([20, 50, 100, 500, 1000])
    ax2.set_yticklabels(["20", "50", "100", "500", "1000"])
    ax2.minorticks_off()

    fig.tight_layout()
    save(fig, "fig17_power_and_fom.png")

    print(f"\n    integrated jitter {j*1e15:.1f} fs at "
          f"{p.p_total_measured*1e3:.1f} mW")
    print(f"    FoM {noise.fom(j, p.p_total_measured):.1f} dB")
    print(f"    FDVPD power {pd_share*p.p_total_measured*1e3:.2f} mW "
          f"({pd_share*100:.0f} % of the total)")


if __name__ == "__main__":
    main()
