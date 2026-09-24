"""Print the design point and its noise budget -- no plotting.

Useful as a first check that an edited parameter set still hangs together, and
as the text version of what the figures show.

    python scripts/design_summary.py
"""

from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from fdvadpll import (default_design, fractional_design, measured_fit_design,
                      noise)

F = np.logspace(3, math.log10(40e6), 4000)


def _jitter(design, frac_bit=None) -> float:
    return noise.integrated_jitter(
        F, noise.output_phase_noise(design, F, frac_bit), design.f_ckv)


def main() -> None:
    d = measured_fit_design()

    print("\nDesign point")
    print("-" * 64)
    print(d.summary())

    print("\nIn-band FDVPD noise budget")
    print("-" * 64)
    print("  integer-N channel:")
    print(noise.pd_noise_floor(d, None))
    print("\n  fractional-N channel (bit 5):")
    print(noise.pd_noise_floor(d, 5))

    print("\nOutput performance")
    print("-" * 64)
    rows = [
        ("pure physics, integer-N", default_design(), None, None),
        ("pure physics, fractional", fractional_design(5), 5, None),
        ("fitted, integer-N", d, None, 82.0),
        ("fitted, fractional (bit 5)", fractional_design(5, d), 5, 101.0),
    ]
    print(f"  {'case':<28} {'jitter':>9}  {'FoM':>9}   measured")
    for label, design, bit, measured in rows:
        j = _jitter(design, bit)
        meas = f"{measured:.0f} fs" if measured else "-"
        print(f"  {label:<28} {j*1e15:7.1f} fs  "
              f"{noise.fom(j, design.power.p_total_measured):7.1f} dB   {meas}")

    print("\nWhy the fit is needed")
    print("-" * 64)
    gap = (10 * math.log10(np.interp(110e3, F,
                                     noise.output_phase_noise(d, F, None)))
           - 10 * math.log10(np.interp(110e3, F,
                                       noise.output_phase_noise(default_design(),
                                                                F, None))))
    print(f"  Sec. V attributes the 1-100 kHz shortfall to ramp-generator")
    print(f"  flicker noise, which the published equations omit.  Adding a")
    print(f"  {d.fdvpd.flicker_corner/1e3:.0f} kHz corner raises L(110 kHz) by "
          f"{gap:.1f} dB and lands on the")
    print(f"  measured -116 dBc/Hz.  No other term is changed.")
    print()


if __name__ == "__main__":
    main()
