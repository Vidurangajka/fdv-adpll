"""MASH sigma-delta modulators.

Used twice in this project:

* dithering the DCO tracking bank so that the effective frequency resolution is
  finer than one ``K_DCO`` LSB (ADPLL extension, Sec. IV of the extension);
* as the *reference* architecture the paper argues against -- a MASH-dithered
  multi-modulus divider, whose shaped quantisation noise is what the FDVPD
  avoids by subtracting the fractional word in the voltage domain (Sec. II-A).
  ``scripts/compare_architectures.py`` uses it for exactly that comparison.
"""

from __future__ import annotations

import numpy as np

__all__ = ["MashSdm", "mash_sequence", "shaped_qn_psd"]


class MashSdm:
    """MASH-1 / 1-1 / 1-1-1 modulator operating on a fractional input.

    ``step(x)`` takes a real value, returns ``round(x)`` plus the shaped
    quantisation error.  With ``n_bits`` the fractional part is first
    quantised to ``2**-n_bits``, modelling a finite-width accumulator.
    """

    def __init__(self, order: int = 1, n_bits: int = 0,
                 rng: np.random.Generator | None = None):
        if order not in (1, 2, 3):
            raise ValueError("order must be 1, 2 or 3")
        self.order = order
        self.n_bits = n_bits
        self.acc = np.zeros(order)
        # Error-feedback delay line.  Always four deep: the recombination below
        # reads _e[1..3] unconditionally, and the unused entries stay zero for
        # the lower orders.
        self._e = np.zeros(4)
        if rng is not None and n_bits > 0:
            # break up idle tones with a half-LSB LSB-dither
            self.acc[0] = rng.uniform(0.0, 1.0)

    def reset(self) -> None:
        self.acc[:] = 0.0
        self._e[:] = 0.0

    def step(self, x: float) -> float:
        integer = np.floor(x)
        frac = x - integer
        if self.n_bits > 0:
            q = 2.0 ** -self.n_bits
            frac = np.floor(frac / q) * q

        carries = np.zeros(self.order)
        u = frac
        for i in range(self.order):
            self.acc[i] += u
            c = np.floor(self.acc[i])
            self.acc[i] -= c
            carries[i] = c
            u = self.acc[i]

        # MASH error-feedback recombination
        y = carries[0]
        if self.order >= 2:
            y += carries[1] - self._e[1]
        if self.order >= 3:
            y += carries[2] - 2.0 * self._e[2] + self._e[3]
        self._e[3] = self._e[2]
        self._e[2] = carries[2] if self.order >= 3 else 0.0
        self._e[1] = carries[1] if self.order >= 2 else 0.0
        return integer + y


def mash_sequence(frac: float, n: int, order: int = 2,
                  rng: np.random.Generator | None = None) -> np.ndarray:
    """Convenience: the integer output sequence for a constant fractional input."""
    m = MashSdm(order=order, rng=rng)
    return np.array([m.step(frac) for _ in range(n)])


def shaped_qn_psd(f, f_ref: float, order: int, t_step: float):
    """Analytic output phase-noise density of a MASH-dithered divider.

    ``L(f) = (2*pi*t_step*f_ref)**2/(12*f_ref) * |2 sin(pi f/f_ref)|**(2*order) / 2``

    with ``t_step`` the time quantisation step (one CKV period for an MMDIV).
    Plotted alongside the FDVPD floor it shows the noise-shaping penalty the
    paper avoids.
    """
    f = np.asarray(f, dtype=float)
    q = (2.0 * np.pi * t_step * f_ref) ** 2 / 12.0 / f_ref
    return 0.5 * q * np.abs(2.0 * np.sin(np.pi * f / f_ref)) ** (2 * order)
