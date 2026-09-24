"""The architectural argument of Sec. II-A, in numbers.

A conventional fractional-N synthesiser dithers a multi-modulus divider with a
MASH modulator: the fractional word is cancelled in the *time* domain, and what
is left is quantisation noise the size of one whole CKV period, shaped upwards
by the modulator.  The FDVPD subtracts the same word in the *voltage* domain
before quantising, so the only quantisation left is the ADC's -- roughly 12 bit
of one PD period.

Two things fall out that are worth seeing side by side:

* deep in band the shaped noise is negligible, so a plot of density alone
  understates the problem;
* a bare type-II loop rolls off at only 20 dB/decade, so a *higher* MASH order
  ends up costing more integrated jitter, not less.  An MMDIV synthesiser needs
  extra filter poles that this architecture does not.

    python scripts/compare_architectures.py
"""

from __future__ import annotations

import math
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from fdvadpll import default_design, noise
from fdvadpll.sdm import shaped_qn_psd
from _style import INK, INK_2, SERIES, save, title, use_style

ORDERS = (1, 2, 3)
FRAC_BIT = 5


def main() -> None:
    use_style()
    d = default_design()
    f = np.logspace(4, math.log10(40e6), 4000)
    h_lp, _ = noise.loop_transfer(d, f)
    shaping = np.abs(h_lp) ** 2

    budget = noise.pd_noise_floor(d, FRAC_BIT)
    fdvpd_q = np.full_like(f, budget.quantisation)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6),
                                  gridspec_kw={"width_ratios": [1.5, 1]})

    # ---- left: densities at the output ---------------------------------
    rows = []
    for order, colour in zip(ORDERS, [SERIES[1], SERIES[3], SERIES[4]]):
        raw = shaped_qn_psd(f, d.f_ref, order, d.t_ckv)
        ax.semilogx(f, 10 * np.log10(raw * shaping), color=colour,
                    label=f"MASH-{order} divider")
        j = noise.integrated_jitter(f, raw * shaping, d.f_ckv)
        rows.append((f"MASH-{order} MMDIV", j * 1e15, colour))

    ax.semilogx(f, 10 * np.log10(fdvpd_q * shaping), color=SERIES[0],
                linewidth=2.4, label="FDVPD quantisation")
    j_fdvpd = noise.integrated_jitter(f, fdvpd_q * shaping, d.f_ckv)
    rows.append(("FDVPD", j_fdvpd * 1e15, SERIES[0]))

    ax.semilogx(f, 10 * np.log10(noise.output_phase_noise(d, f, FRAC_BIT)),
                color=INK, linewidth=1.4, linestyle=(0, (4, 3)),
                label="full FDVPD design (all sources)")

    title(ax, "Phase quantisation noise at the output",
          "both after the same type-II loop filter")
    ax.set_xlabel("offset frequency [Hz]")
    ax.set_ylabel("L(f) [dBc/Hz]")
    ax.set_ylim(-190, -100)
    # lower left: the only corner the rising MASH traces never reach
    ax.legend(loc="lower left")

    # ---- right: what it costs in jitter ---------------------------------
    rows.sort(key=lambda r: r[1])
    names = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    colours = [r[2] for r in rows]
    y = np.arange(len(rows))
    ax2.barh(y, vals, color=colours, height=0.6)
    ax2.set_yticks(y)
    ax2.set_yticklabels(names)
    ax2.invert_yaxis()
    ax2.grid(axis="y", visible=False)
    for yi, v in zip(y, vals):
        ax2.annotate(f"{v:.0f} fs", xy=(v, yi), xytext=(6, 0),
                     textcoords="offset points", va="center", fontsize=8.5,
                     color=INK_2)
    title(ax2, "Jitter each one contributes", "rms, 10 kHz - 40 MHz")
    ax2.set_xlabel("jitter contribution [fs]")
    ax2.set_xlim(0, max(vals) * 1.25)

    fig.tight_layout()
    save(fig, "compare_architectures.png")

    print()
    for name, v, _ in rows:
        print(f"    {name:16s} {v:7.1f} fs")
    print(f"\n    the paper's whole measured budget is 82 fs (integer-N), "
          f"101 fs (fractional-N)")


if __name__ == "__main__":
    main()
