"""Behavioural models of the FDVPD blocks (Fig. 6).

    10 b segmented current DAC  ->  differential dv/dt ramp  ->  7 b SAR ADC

Each block carries the imperfections that Sec. III-E lists as the linearity
limiters of the phase detector, so the fractional-spur mechanisms emerge from
the simulation instead of being added by hand:

    1) I-DAC static INL/DNL and code-dependent settling
    2) ramp gain error (I_R mismatch) and 2nd/3rd-order ramp non-linearity
    3) switching charge injection at S1 / S2
    4) CDAC mismatch inside the SAR ADC

Timing convention (Fig. 7(b))
-----------------------------
Let ``Phi_R[k] = k * FCW_pd`` be the ideal accumulated phase of CKVd, in CKVd
cycles, at reference edge ``k``.  The gated edge CKVdg is the first CKVd edge
at or after REF, so its ideal delay is

    dt_ideal[k] = T_frac[k] + T_const,
    T_frac[k]   = (ceil(Phi_R[k]) - Phi_R[k]) * T_PD

which is equation (14): for a non-zero fractional part ``ceil(x) - x`` is
``1 - {x}``.  ``T_const = V_OS/SR`` is the fixed margin that keeps REF ahead of
the expected edge.  In an integer-N channel ``T_frac`` collapses to zero and
the ramp window is just ``T_const`` -- the "fixed ON-time of 100 ps" that the
paper uses to model the integer-N ``I_R`` contribution in Fig. 12.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .params import KT, DesignParams, FdvpdParams

__all__ = ["CurrentDAC", "RampGenerator", "SarAdc", "Fdvpd", "FdvpdSample"]


# --------------------------------------------------------------------------
class CurrentDAC:
    """10 b segmented I-DAC: 4 b thermometer MSB + 6 b binary LSB (Sec. III-C).

    Unit-element mismatch uses the usual ``sigma ~ sqrt(W)`` law, so a
    thermometer element of weight 64 LSB has 8x the sigma of the 1 LSB element.
    The whole 2**bits transfer characteristic is built once, which reduces the
    per-cycle cost to an array lookup.
    """

    def __init__(self, params: FdvpdParams, rng: np.random.Generator):
        self.bits = params.dac_bits
        self.thermo_bits = params.dac_msb_thermo_bits
        self.bin_bits = self.bits - self.thermo_bits
        self.n_codes = 2 ** self.bits

        seg = 2 ** self.bin_bits                       # weight of one MSB unit
        n_thermo = 2 ** self.thermo_bits - 1
        w_thermo = np.full(n_thermo, float(seg))
        w_bin = 2.0 ** np.arange(self.bin_bits)        # 1, 2, 4, ..., seg/2

        sig = params.dac_sigma_lsb
        if sig > 0.0:
            w_thermo = w_thermo + rng.standard_normal(n_thermo) * sig * np.sqrt(w_thermo)
            w_bin = w_bin + rng.standard_normal(self.bin_bits) * sig * np.sqrt(w_bin)
        self.w_thermo, self.w_bin = w_thermo, w_bin

        codes = np.arange(self.n_codes)
        msb = codes >> self.bin_bits
        lsb = codes & (seg - 1)
        cum = np.concatenate(([0.0], np.cumsum(w_thermo)))   # thermometer sum
        bits_lsb = (lsb[:, None] >> np.arange(self.bin_bits)) & 1
        self.levels = cum[msb] + bits_lsb @ w_bin            # LSB units
        self.levels *= (1.0 + params.dac_gain_err)
        # fraction of one PD period represented by each code
        self.frac = self.levels / self.n_codes

    def _endpoint_fit(self):
        x = np.arange(self.n_codes, dtype=float)
        gain = (self.levels[-1] - self.levels[0]) / (self.n_codes - 1)
        return self.levels[0] + gain * x

    @property
    def inl(self) -> np.ndarray:
        """Integral non-linearity [LSB], endpoint fit."""
        return self.levels - self._endpoint_fit()

    @property
    def dnl(self) -> np.ndarray:
        """Differential non-linearity [LSB], endpoint fit."""
        gain = (self.levels[-1] - self.levels[0]) / (self.n_codes - 1)
        return np.diff(self.levels) / gain - 1.0


# --------------------------------------------------------------------------
class RampGenerator:
    """Differential dv/dt ramp, ``SR = 2*I_R/C_SAR`` (Sec. III-B).

    * ``ramp_mismatch`` -- mismatch between the top and bottom ``I_R``.  Because
      CKVdg always lags REF, both sources carry the *same* lag information, so
      the mismatch appears only as a static gain error rather than as the
      lead/lag asymmetry that produces charge-pump spurs (Sec. III-A).
    * ``ramp_nl2/nl3`` -- residual channel-length modulation, applied as a
      memoryless distortion of the ramp excursion.  A constant offset and a
      linear gain error are absorbed by the loop, so only 2nd order and above
      generate spurs.
    """

    def __init__(self, params: FdvpdParams):
        self.p = params
        self.sr = params.slew_rate * (1.0 + params.ramp_mismatch)
        s_i = 4.0 * KT * params.gamma * 2.0 * params.gm_ramp   # equation (2)
        self._i_noise_k = math.sqrt(s_i / 2.0) / params.c_sar
        self.ktc_sigma = (math.sqrt(2.0 * KT / params.c_sar)   # equation (9)
                          if params.ktc_enabled else 0.0)

    def excursion(self, dt):
        """Differential voltage swept in time ``dt`` [s]."""
        e = self.sr * np.asarray(dt, dtype=float)
        return e * (1.0 + self.p.ramp_nl2 * e + self.p.ramp_nl3 * e * e)

    def current_noise_sigma(self, t_on):
        """RMS sampled voltage noise from ``I_R`` integrated over ``t_on``.

        White current noise of one-sided PSD ``S`` integrated over a window
        ``t_on`` has charge variance ``S*t_on/2``; dividing by ``C_SAR`` gives
        the sampled voltage.  Substituting into ``sigma_t = sigma_v/SR`` and
        converting to phase reproduces equation (7) exactly.
        """
        return self._i_noise_k * np.sqrt(np.maximum(np.asarray(t_on), 0.0))


# --------------------------------------------------------------------------
class SarAdc:
    """7 b top-plate-sampling monotonic-switching self-timed SAR (Sec. III-D).

    With static capacitor mismatch a binary-search SAR is still a memoryless
    quantiser, so its behaviour is fully captured by the 2**bits decision
    levels, built once and then applied with ``searchsorted``.
    """

    def __init__(self, params: FdvpdParams, rng: np.random.Generator,
                 gain_trim: float = 1.0):
        self.bits = params.adc_bits
        self.n_codes = 2 ** self.bits
        self.lsb = params.adc_lsb * gain_trim
        self.v_noise = params.comparator_noise

        w = 2.0 ** np.arange(self.bits)
        if params.adc_sigma_cap > 0.0:
            w = w * (1.0 + rng.standard_normal(self.bits)
                     * params.adc_sigma_cap / np.sqrt(w))
        codes = np.arange(self.n_codes)
        bits_m = (codes[:, None] >> np.arange(self.bits)) & 1
        self.levels = ((bits_m @ w) - self.n_codes / 2.0) * self.lsb
        if not np.all(np.diff(self.levels) > 0):
            raise ValueError("CDAC mismatch made the SAR non-monotonic; "
                             "reduce adc_sigma_cap")

    @property
    def full_scale(self) -> float:
        return self.n_codes * self.lsb

    def convert(self, v_in, noise=None):
        """Digitise ``v_in`` [V]; ``noise`` is a pre-drawn unit normal (or array).

        Returns ``(signed_code, saturated)``.  ``saturated`` marks the
        bang-bang regime the loop sits in before lock (Sec. III-A, step 5).
        """
        v = np.asarray(v_in, dtype=float)
        if noise is not None and self.v_noise > 0.0:
            v = v + noise * self.v_noise
        idx = np.searchsorted(self.levels, v, side="right") - 1
        sat = (idx < 0) | (idx >= self.n_codes - 1)
        return np.clip(idx, 0, self.n_codes - 1) - self.n_codes // 2, sat


# --------------------------------------------------------------------------
@dataclass
class FdvpdSample:
    """One phase-detector conversion."""

    code: int           # signed SAR output
    error_time: float   # code translated back to a time error [s]
    v_sampled: float    # differential voltage at the ADC input [V]
    t_on: float         # ramp window [s]
    dac_code: int
    saturated: bool


class Fdvpd:
    """Fully differential voltage phase detector, steps 1-5 of Sec. III-A.

      step 1  encode   -- DAC drives ``V_DM = SR*T_frac`` onto C_SAR
      step 2  transfer -- REF starts the ramp, C_SAR discharges at ``-SR``
      step 3  sample   -- CKVdg stops the ramp and samples the residue
      step 4  extend   -- ``C_exten`` adds a fixed offset, i.e. a fixed phase
                          offset that the loop regulates away (Sec. III-D)
      step 5  convert  -- self-timed SAR digitises the residue

    Call :meth:`prime` once with the number of reference cycles to pre-draw the
    noise; ``step`` then costs a handful of scalar operations.
    """

    def __init__(self, design: DesignParams, rng: np.random.Generator,
                 gain_trim: float = 1.0):
        self.design = design
        self.p = design.fdvpd
        self.rng = rng
        self.dac = CurrentDAC(self.p, rng)
        self.ramp = RampGenerator(self.p)
        self.adc = SarAdc(self.p, rng, gain_trim)
        self.t_pd = design.t_pd
        self._v_dac_prev = 0.0
        self._n = 0
        self._k = 0
        self._z_ramp = self._z_ktc = self._z_cmp = None
        self._v_flicker = None
        # gain from a sampled voltage back to a time error; the calibration
        # loop in calib.py trims this at run time.
        self.time_per_lsb = self.p.adc_lsb / self.p.slew_rate

    # ------------------------------------------------------------------
    def prime(self, n: int, mean_t_on: float | None = None) -> None:
        """Pre-draw ``n`` samples of every phase-detector noise source.

        The 1/f term needs the *white* voltage variance to set its level, so
        ``mean_t_on`` (the average ramp window, equation (1)) is required when
        ``flicker_corner`` is non-zero.
        """
        from .dsp import one_over_f

        self._n = int(n)
        self._k = 0
        self._z_ramp = self.rng.standard_normal(self._n)
        self._z_ktc = self.rng.standard_normal(self._n)
        self._z_cmp = self.rng.standard_normal(self._n)

        fc = self.p.flicker_corner
        if fc > 0.0:
            t_on = self.t_pd / 2.0 if mean_t_on is None else mean_t_on
            v_white_sq = (self.ramp.current_noise_sigma(t_on) ** 2
                          + self.ramp.ktc_sigma ** 2
                          + self.p.comparator_noise ** 2)
            # A white voltage sequence of variance v^2 at f_REF has one-sided
            # PSD 2v^2/f_REF; the 1/f part must cross it at f = fc, so its PSD
            # is (2 v^2 / f_REF) * fc / f.
            f_ref = self.design.f_ref
            amp = 2.0 * float(v_white_sq) * fc / f_ref
            self._v_flicker = one_over_f(self._n, f_ref, amp, self.rng)
        else:
            self._v_flicker = None

    def _z(self, arr):
        if self._k < self._n:
            return arr[self._k]
        return self.rng.standard_normal()

    # ------------------------------------------------------------------
    def dac_code(self, t_frac_norm: float) -> int:
        """Code for ``T_frac`` expressed as a fraction of one PD period."""
        code = int(round(t_frac_norm * self.dac.n_codes))
        return min(max(code, 0), self.dac.n_codes - 1)

    def step(self, dt: float, t_frac_norm: float,
             charge_inj: float = 0.0) -> FdvpdSample:
        """One conversion.

        ``dt``           measured REF -> CKVdg delay [s]
        ``t_frac_norm``  predicted ``T_frac / T_PD``.  The ``V_OS`` margin is
                         already folded into it by the reference accumulator
                         (see :mod:`fdvadpll.pll`), so the lock point is simply
                         ``dt = t_frac_norm * T_PD``.
        """
        p = self.p
        code = self.dac_code(t_frac_norm)
        v_dm = self.ramp.sr * self.dac.frac[code] * self.t_pd

        if p.dac_settle_tau > 0.0:      # code-dependent settling (Sec. III-E1)
            t_avail = 1.0 / self.design.f_ref - dt
            k = math.exp(-max(t_avail, 0.0) / p.dac_settle_tau)
            v_dm = v_dm + (self._v_dac_prev - v_dm) * k
        self._v_dac_prev = v_dm

        v = v_dm - self.ramp.excursion(dt)                    # steps 2 + 3
        v = float(v)
        v += self._z(self._z_ramp) * self.ramp.current_noise_sigma(dt)
        v += self._z(self._z_ktc) * self.ramp.ktc_sigma
        if self._v_flicker is not None:
            v += self._z(self._v_flicker)
        v += charge_inj

        c, sat = self.adc.convert(v, self._z(self._z_cmp))    # steps 4 + 5
        self._k += 1
        c = int(c)
        return FdvpdSample(c, c * self.time_per_lsb, v, float(dt), code,
                           bool(sat))
