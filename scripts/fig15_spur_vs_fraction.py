"""Fig. 15 -- fractional spur versus the asserted fractional bit.

Sweeps the near-integer channels FCW_pd = 21 - 2^-Nq and measures the spur at
the fractional frequency.  Three phase detectors are compared: the ideal one,
one with a 0.3 % DAC-to-ramp gain mismatch, and the same mismatch with the
background gain calibration running.

Deep fractional bits move the spur inside the loop bandwidth, where the loop no
longer attenuates it -- which is why the ideal trace rises towards the right
rather than staying flat.

On the size of the mismatch: a *positive* gain error makes the DAC ask for more
than one PD period at the top of the sawtooth, so the top of the fractional
range is simply unreachable and the detector rails once per sawtooth.  Below
about 0.5 % the loop rides through that; at 1 % it loses lock outright in the
deep channels, and then there is nothing for a background calibration to adapt
on.  0.3 % is the interesting regime -- large enough to dominate the spur,
small enough that the loop still works.

    python scripts/fig15_spur_vs_fraction.py
"""

from __future__ import annotations

import sys
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from fdvadpll import FdvPll, fractional_design, measured_fit_design
from _style import INK_2, SERIES, save, title, use_style

BITS = list(range(2, 12))
N_CYCLES = 1 << 16
SEED = 5

GAIN_ERR = 0.003

#: (legend label, short direct label, PD overrides, run overrides, colour)
CASES = [
    ("ideal I-DAC", "ideal", {}, {}, SERIES[0]),
    (f"{GAIN_ERR*100:.1f} % gain mismatch", "mismatch",
     {"dac_gain_err": GAIN_ERR}, {}, SERIES[1]),
    ("mismatch + gain cal.", "calibrated", {"dac_gain_err": GAIN_ERR},
     {"calibrate": ("gain",), "cal_mu": {"gain": 5e-3}}, SERIES[2]),
]


def main() -> None:
    use_style()
    base = measured_fit_design()
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    for label, short, pd_kw, run_kw, colour in CASES:
        spurs, jitters = [], []
        for bit in BITS:
            d = fractional_design(bit, base)
            if pd_kw:
                d = replace(d, fdvpd=replace(d.fdvpd, **pd_kw))
            res = FdvPll(d, seed=SEED, **run_kw).run(N_CYCLES)
            spurs.append(res.fractional_spur()[1])
            jitters.append(res.jitter() * 1e15)
            print(f"    {label:22s} bit {bit:2d}: "
                  f"{spurs[-1]:7.1f} dBc, {jitters[-1]:6.1f} fs")
        ax.plot(BITS, spurs, color=colour, marker="o", label=label)
        ax2.plot(BITS, jitters, color=colour, marker="o", label=label)
        ax.annotate(short, xy=(BITS[-1], spurs[-1]), xytext=(7, 0),
                    textcoords="offset points", fontsize=8.5, color=INK_2,
                    va="center")

    title(ax, "Fractional spur",
          f"FCW_pd = 21 - 2^-Nq,  {N_CYCLES} reference cycles per point")
    ax.set_xlabel("asserted fractional bit  Nq")
    ax.set_ylabel("spur at f_frac [dBc]")
    ax.set_xticks(BITS)
    ax.set_xlim(BITS[0] - 0.4, BITS[-1] + 1.6)   # room for the direct labels
    ax.legend(loc="lower right")

    title(ax2, "Integrated jitter",
          "rms, 10 kHz - 40 MHz; log scale -- channels that lose lock run to "
          "nanoseconds")
    ax2.set_xlabel("asserted fractional bit  Nq")
    ax2.set_ylabel("jitter [fs]")
    ax2.set_yscale("log")
    ax2.set_xticks(BITS)
    ax2.legend(loc="upper left")

    fig.tight_layout()
    save(fig, "fig15_spur_vs_fraction.png")


if __name__ == "__main__":
    main()
