"""Fig. 11 -- measured output phase noise, integer-N and fractional-N.

Plots the event-driven simulation against the independent s-domain model of
Sec. II-B.  The two share nothing but the design parameters, so their agreement
is the main evidence that either one is right.

The published anchor points (Sec. V and the Fig. 11 annotations) are drawn as
markers.  The full measured trace is not reproduced here -- the paper does not
publish it as data -- so only the numbers it actually quotes are shown.

    python scripts/fig11_phase_noise.py
"""

from __future__ import annotations

import math
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from fdvadpll import FdvPll, fractional_design, measured_fit_design, noise
from _style import SERIES, annotate, ordinal, save, title, use_style

N_CYCLES = 1 << 17
FRAC_BIT = 5

# Published anchor points (Wu et al., JSSC 2021)
MEASURED = {
    "integer": {"jitter_fs": 82.0, "at_110k_dbc": -116.0},
    "fractional": {"jitter_fs": 101.0},
}


def main() -> None:
    use_style()
    f_grid = np.logspace(3, math.log10(40e6), 4000)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)

    for ax, (label, design, frac_bit) in zip(axes, [
        ("Integer-N  (FCW = 21 x 2)", measured_fit_design(), None),
        (f"Fractional-N  (FCW = 20.9688 x 2)",
         fractional_design(FRAC_BIT, measured_fit_design()), FRAC_BIT),
    ]):
        print(f"  simulating {label} ...")
        res = FdvPll(design, seed=11).run(N_CYCLES)
        f_sim, L_sim = res.phase_noise()
        L_ana = noise.output_phase_noise(design, f_grid, frac_bit)

        ax.semilogx(f_sim, 10 * np.log10(L_sim), color=SERIES[0], linewidth=1.1,
                    alpha=0.75, label="simulated (event driven)")
        ax.semilogx(f_grid, 10 * np.log10(L_ana), color=SERIES[1],
                    label="analytic (Sec. II-B)")

        j_sim = res.jitter() * 1e15
        j_ana = noise.integrated_jitter(f_grid, L_ana, design.f_ckv) * 1e15
        key = "integer" if frac_bit is None else "fractional"
        j_meas = MEASURED[key]["jitter_fs"]

        if "at_110k_dbc" in MEASURED[key]:
            ax.plot([110e3], [MEASURED[key]["at_110k_dbc"]], marker="D",
                    color=SERIES[2], markeredgecolor="#fcfcfb",
                    markeredgewidth=1.5, linestyle="none", zorder=6,
                    label="measured (published)")
            annotate(ax, 110e3, MEASURED[key]["at_110k_dbc"],
                     f"{MEASURED[key]['at_110k_dbc']:.0f} dBc/Hz  ",
                     SERIES[2], dx=-6, dy=10, ha="right", va="bottom")

        title(ax, label,
              f"sim {j_sim:.0f} fs  |  analytic {j_ana:.0f} fs  |  "
              f"measured {j_meas:.0f} fs   (rms, 10 kHz - 40 MHz)")
        ax.set_xlabel("offset frequency [Hz]")
        ax.set_xlim(1e4, 4e7)
        ax.set_ylim(-160, -95)
        ax.legend(loc="upper right")

    axes[0].set_ylabel("L(f) [dBc/Hz]")
    fig.tight_layout()
    save(fig, "fig11_phase_noise.png")


if __name__ == "__main__":
    main()
