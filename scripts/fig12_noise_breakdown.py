"""Fig. 12 -- where the output phase noise comes from.

Each contribution is the source multiplied by the loop transfer function it
actually sees: the reference and the four phase-detector terms through the
low-pass, the oscillator through the high-pass.  The free-running DCO trace is
drawn dashed so the loop's rejection in band is visible.

    python scripts/fig12_noise_breakdown.py
"""

from __future__ import annotations

import math
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from fdvadpll import fractional_design, measured_fit_design, noise
from _style import INK, MUTED, SERIES, save, title, use_style

FRAC_BIT = 5

#: fixed slot per contribution -- colour follows the entity, not its rank
COLOURS = {
    "reference": SERIES[0],
    "ramp current": SERIES[1],
    "kT/C": SERIES[2],
    "comparator": SERIES[3],
    "quantisation": SERIES[4],
    "DCO": SERIES[5],
}


def panel(ax, design, frac_bit, label) -> None:
    f = np.logspace(3, math.log10(40e6), 4000)
    total, parts = noise.output_phase_noise(design, f, frac_bit,
                                            breakdown=True)

    for name, colour in COLOURS.items():
        ax.semilogx(f, 10 * np.log10(parts[name]), color=colour, linewidth=1.6,
                    label=name)
    ax.semilogx(f, 10 * np.log10(parts["DCO (free running)"]), color=SERIES[5],
                linewidth=1.2, linestyle=(0, (4, 3)), alpha=0.65,
                label="DCO (free running)")
    ax.semilogx(f, 10 * np.log10(total), color=INK, linewidth=2.4,
                label="total")

    ax.axvline(design.loop.bandwidth, color=MUTED, linewidth=0.8,
               linestyle=(0, (2, 3)))
    ax.text(design.loop.bandwidth * 1.12, -157,
            f"loop BW\n{design.loop.bandwidth/1e3:.0f} kHz", fontsize=8,
            color=MUTED, va="bottom")

    jitter = noise.integrated_jitter(f, total, design.f_ckv)
    title(ax, label, f"{jitter*1e15:.0f} fs rms  |  FoM "
                     f"{noise.fom(jitter, design.power.p_total_measured):.1f} dB")
    ax.set_xlabel("offset frequency [Hz]")
    ax.set_xlim(1e4, 4e7)
    ax.set_ylim(-175, -95)


def main() -> None:
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    panel(axes[0], measured_fit_design(), None, "Integer-N")
    panel(axes[1], fractional_design(FRAC_BIT, measured_fit_design()),
          FRAC_BIT, "Fractional-N  (FCW = 20.9688 x 2)")
    axes[0].set_ylabel("L(f) [dBc/Hz]")
    # lower left is the only region no trace passes through
    axes[1].legend(loc="lower left", ncol=2)
    fig.tight_layout()
    save(fig, "fig12_noise_breakdown.png")

    # the in-band budget, as a table
    print("\n  In-band FDVPD noise budget (fractional, bit 5):")
    print(noise.pd_noise_floor(fractional_design(FRAC_BIT,
                                                 measured_fit_design()),
                               FRAC_BIT))


if __name__ == "__main__":
    main()
