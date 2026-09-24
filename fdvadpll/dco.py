"""Digitally controlled oscillator.

Two tuning models are provided:

``ideal``
    ``f = f_center + NTW * f_REF``.  The normalised tuning word feeds the
    oscillator directly, i.e. ``K_DCO`` is perfectly known.  This is the model
    used to reproduce the paper's measurements.

``banks``
    The ADPLL extension: three physical tuning banks (PVT / acquisition /
    tracking) with integer oscillator tuning words, the tracking bank dithered
    by a MASH sigma-delta to gain sub-LSB resolution.  The loop converts its
    normalised tuning word with an *estimated* ``K_DCO``, so a ``kdco_err``
    directly scales the loop gain -- which is what the ``KdcoCalibration`` in
    ``calib.py`` removes.

Phase noise is injected as a per-edge period deviation, which keeps the
oscillator's own noise independent of the loop and lets the loop act on it
exactly as it does in hardware.
"""

from __future__ import annotations

import math

import numpy as np

from .dsp import dco_period_jitter
from .params import DesignParams
from .sdm import MashSdm

__all__ = ["Dco"]


class Dco:
    """Harmonic-shaped LC DCO (Fig. 9) running at the divided (CKVd) rate.

    The simulator advances CKVd edges, so the period jitter is referred to
    CKVd: summing ``fb_div`` CKV periods multiplies the period variance by
    ``fb_div`` and, at the ``fb_div``-times lower carrier, lowers the phase
    noise by ``20*log10(fb_div)`` -- the standard divider result.
    """

    def __init__(self, design: DesignParams, rng: np.random.Generator,
                 n_edges: int, mode: str = "ideal", noise: bool = True):
        if mode not in ("ideal", "banks"):
            raise ValueError("mode must be 'ideal' or 'banks'")
        self.design = design
        self.p = design.dco
        self.mode = mode
        self.div = design.fb_div
        self.f_center = self.p.f_center

        self.f_edge_nom = design.f_pd            # CKVd rate
        sigma_ckv = self.p.period_jitter_white()
        self.sigma_edge = sigma_ckv * math.sqrt(self.div)
        if noise:
            self.dt_noise = dco_period_jitter(
                n_edges + 8, self.f_edge_nom, self.sigma_edge,
                self.p.f_corner_flicker, rng)
        else:
            self.dt_noise = np.zeros(n_edges + 8)
        self._edge = 0

        # --- tuning state -------------------------------------------------
        self.f_dco = self.p.f_center
        self.ntw = 0.0
        self.otw_pvt = 0
        self.otw_acq = 0
        self.otw_trk = 0.0
        self.kdco_hat = self.p.kdco * (1.0 + self.p.kdco_err)
        self._dither = MashSdm(order=1, n_bits=self.p.trk_dither_bits, rng=rng)
        self._trk_int = 0

    # ------------------------------------------------------------------
    @property
    def kdco_true(self) -> float:
        return self.p.kdco

    @property
    def coarse_hz(self) -> float:
        """Frequency the PVT and acquisition banks contribute [Hz].

        Their gains are fixed design constants, so ``KdcoCalibration`` can
        subtract this exactly and regress only what the tracking bank did.
        """
        return (self.otw_pvt * self.p.kdco_pvt
                + self.otw_acq * self.p.kdco_acq)

    def set_normalised_tuning(self, ntw: float) -> None:
        """Apply a normalised tuning word (``delta_f = NTW * f_REF`` if ideal)."""
        self.ntw = ntw
        if self.mode == "ideal":
            self.f_dco = self.f_center + ntw * self.design.f_ref
            return
        # --- banks: convert with the *estimated* K_DCO --------------------
        want = ntw * self.design.f_ref / self.kdco_hat      # in tracking LSBs
        lo, hi = -2 ** (self.p.trk_bits - 1), 2 ** (self.p.trk_bits - 1) - 1
        acq_per_trk = self.p.kdco_acq / self.p.kdco

        # ``want`` is the *total* offset, recomputed from scratch every cycle,
        # so the coarse bank's share comes off every cycle -- not only on the
        # cycle that moved it.  The coarse bank is then re-centred whenever the
        # residue no longer fits the tracking bank, which is what lets it come
        # back down again as well as go up; leaving it where acquisition
        # happened to put it would strand the oscillator megahertz off target.
        residue = want - self.otw_acq * acq_per_trk
        if residue > hi or residue < lo:
            self.otw_acq = int(np.clip(
                round(want / acq_per_trk),
                -2 ** (self.p.acq_bits - 1), 2 ** (self.p.acq_bits - 1) - 1))
            residue = want - self.otw_acq * acq_per_trk
        self.otw_trk = float(np.clip(residue, lo, hi))
        coded = self._dither.step(self.otw_trk)              # integer + dither
        self.f_dco = (self.f_center
                      + self.otw_pvt * self.p.kdco_pvt
                      + self.otw_acq * self.p.kdco_acq
                      + coded * self.p.kdco)

    def set_frequency(self, f: float) -> None:
        """Force the instantaneous frequency (used to seed the initial state)."""
        self.f_dco = f

    # ------------------------------------------------------------------
    def periods(self, n: int) -> np.ndarray:
        """Next ``n`` CKVd periods [s], including the oscillator's own noise."""
        t0 = self.div / self.f_dco
        k = self._edge
        if k + n > len(self.dt_noise):      # grow *before* slicing, or the
            extra = max(k + n - len(self.dt_noise), n * 4, 1024)   # slice comes
            self.dt_noise = np.concatenate(  # back short and the caller's edge
                [self.dt_noise, np.zeros(extra)])                  # times drift
        out = t0 + self.dt_noise[k:k + n]
        self._edge += n
        return out

    def free_running_pn(self, f):
        """Free-running L(f) [linear] at the *output* carrier ``f_CKV``."""
        f = np.asarray(f, dtype=float)
        lin = 10.0 ** (self.p.pn_1mhz / 10.0)
        white = lin * (1e6 / f) ** 2
        return white * (1.0 + self.p.f_corner_flicker / f)
