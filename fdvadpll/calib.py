"""Background digital calibration loops -- the ADPLL extension.

Sec. III of the paper notes that the two dominant static imperfections of the
FDVPD ("a mismatch between the I_R sources only results in static deviations
from the expected SR, which can easily be calibrated", and the I-DAC INL that
"motivates the adoption of a segmented I-DAC") are *calibratable*, but does not
give the loops.  This module supplies them:

``PdGainCalibration``
    Corrects the gain mismatch between the voltage the I-DAC produces for
    ``T_frac`` and the Volts-per-second of the ramp.  This mismatch makes the
    locked edge land at ``(SR_dac/SR_ramp)*T_frac`` instead of ``T_frac``, i.e.
    it turns the fractional accumulator's sawtooth into an output phase
    modulation -- the dominant fractional-spur mechanism.

``DacInlCalibration``
    A 15-entry look-up table indexed by the thermometer MSB segment, adapted by
    a sign-LMS.  It removes the segment-boundary INL that shows up as harmonics
    of the fractional frequency.

``KdcoCalibration``
    Tracks ``K_DCO`` from the observed frequency step per unit tuning word, so
    the normalised loop coefficients keep their designed bandwidth over PVT.

Why this converges
------------------
In lock the loop drives the *measured* phase error towards zero, so a naive
correlation would starve.  It does not, because the loop only suppresses the
error inside its bandwidth: the residue seen by the ADC is ``H_hp(s)`` times
the injected error, and the fractional sawtooth has harmonics above the loop
bandwidth.  Convergence is therefore fastest when the calibration runs at a
*reduced* loop bandwidth -- the opposite of the usual intuition, and worth
remembering when porting this to silicon.
"""

from __future__ import annotations

import numpy as np

__all__ = ["PdGainCalibration", "DacInlCalibration", "KdcoCalibration"]


class PdGainCalibration:
    """LMS trim of the DAC-to-ramp gain ratio.

    The correction multiplies the fractional word before it reaches the DAC:

        T_frac_applied = g_hat * T_frac_ideal

    ``update`` is called once per reference cycle with the measured phase error
    and the fractional word that produced it -- and only when the ADC is in
    range, since a railed code carries no gradient.

    Two guard rails, because the gradient is not always trustworthy.  The
    regressor is the fractional sawtooth; once that sawtooth falls *inside* the
    loop bandwidth, what reaches the detector is the loop's high-pass-shaped
    version of it, which is phase-rotated.  Correlating a rotated signal
    against the unrotated one scales the gradient by the cosine of that
    rotation -- and past 90 degrees the sign inverts, so in the deepest
    channels the loop can adapt confidently in the wrong direction.  (This is
    the same effect the module docstring describes from the other side: the
    calibration converges best at a *reduced* loop bandwidth.)

    ``clamp`` therefore has to be a physical limit, not a generous one.  The
    detector can absorb a lock-point error of ``adc_full_scale/SR`` out of one
    PD period -- about 3 % for the default design -- so a clamp wider than that
    lets a wrong estimate rail the ADC, which stops the updates, which removes
    any chance of recovery.

    ``leak`` pulls the estimate gently back towards unity and, unlike the
    correlation term, keeps running while the detector is railed (call
    :meth:`relax` on those cycles).  That is what makes a bad excursion
    recoverable rather than permanent.  It is small enough not to measurably
    bias a channel that does have a usable gradient.
    """

    def __init__(self, mu: float = 2e-3, g_init: float = 1.0,
                 sign_sign: bool = False, clamp: float = 0.03,
                 leak: float = 1e-5):
        self.mu = mu
        self.g = g_init
        self.sign_sign = sign_sign
        self.clamp = clamp
        self.leak = leak
        self.history: list[float] = []

    def correct(self, t_frac_norm: float) -> float:
        return float(np.clip(self.g * t_frac_norm, 0.0, 1.0 - 1e-12))

    def update(self, phi_e: float, t_frac_norm: float) -> None:
        x = t_frac_norm - 0.5          # zero-mean regressor
        if self.sign_sign:
            self.g += self.mu * np.sign(phi_e) * np.sign(x)
        else:
            self.g += self.mu * phi_e * x
        self.relax()

    def relax(self) -> None:
        """Apply the leak alone -- for cycles with no usable measurement.

        Called on every railed sample, so an excursion that stopped the
        updates still decays instead of sticking at the clamp forever.
        """
        self.g -= self.leak * (self.g - 1.0)
        self.g = float(np.clip(self.g, 1.0 - self.clamp, 1.0 + self.clamp))
        self.history.append(self.g)


class DacInlCalibration:
    """Sign-LMS look-up table indexed by the I-DAC thermometer segment.

    One entry per MSB segment (15 for a 4 b thermometer array).  The correction
    is applied in the fractional-word domain, in units of one DAC LSB.
    """

    def __init__(self, n_segments: int = 16, mu: float = 2e-4,
                 dac_bits: int = 10, bin_bits: int = 6):
        self.lut = np.zeros(n_segments)
        self.mu = mu
        self.bin_bits = bin_bits
        self.scale = 1.0 / (2 ** dac_bits)     # one DAC LSB, normalised

    def segment(self, dac_code: int) -> int:
        return int(dac_code) >> self.bin_bits

    def correct(self, t_frac_norm: float, dac_code: int) -> float:
        s = self.segment(dac_code)
        return float(np.clip(t_frac_norm + self.lut[s] * self.scale,
                             0.0, 1.0 - 1e-12))

    def update(self, phi_e: float, dac_code: int) -> None:
        s = self.segment(dac_code)
        self.lut[s] += self.mu * np.sign(phi_e)
        self.lut -= self.lut.mean()            # keep it zero-mean (no DC term)


class KdcoCalibration:
    """Running estimate of ``K_DCO`` from observed frequency steps.

    The counter read-out is a frequency meter: the CKVd count accumulated over
    ``window`` reference periods gives

        f_hat = (count[k] - count[k-window]) * f_REF * div / window

    Differencing two consecutive windows removes the unknown centre frequency
    and leaves a frequency *step*, which is regressed on the matching step in
    the (window-averaged) tuning word by a plain LMS:

        K_DCO_hat += mu * (df - K_DCO_hat * d_otw) * d_otw / d_otw^2

    Sizing the window is the whole design problem.  The counter is an integer,
    so one window resolves frequency only to ``f_REF*div/window`` -- 160 MHz
    divided by the window length here.  The tracking bank spans just
    ``2**trk_bits * K_DCO`` = 3.1 MHz end to end, so a short window cannot see
    a tracking-word step at all: every measurement returns the same count, the
    apparent ``df`` is zero, and the LMS obligingly drives ``K_DCO_hat`` to
    zero.  The window must therefore be long enough that

        f_REF * div / window  <<  K_DCO * (expected step in LSBs)

    which for the default design means tens of thousands of reference cycles,
    i.e. a fraction of a millisecond.  That is a property of counter-based
    frequency measurement, not of this loop; a silicon implementation with a
    higher-resolution frequency detector can use a much shorter one.  Steps
    smaller than ``min_step`` tuning LSBs are skipped, since there the
    measurement is all quantisation noise.

    ``coarse_hz`` is the frequency the PVT and acquisition banks contribute.
    It must be supplied, because those banks move during acquisition -- exactly
    when the tracking word is slewing enough to be worth regressing on -- and a
    frequency step they caused would otherwise be credited to ``K_DCO``.  Their
    own gains are fixed design constants, so subtracting them is exact.

    In lock the tuning word barely moves, so there is nothing to learn from;
    essentially all of the useful updates happen while the loop is still
    acquiring.  That is the normal situation for this kind of background loop,
    not a limitation of the implementation.
    """

    def __init__(self, kdco_init: float, mu: float = 0.25,
                 min_step: float = 0.5, window: int = 16384):
        self.kdco = kdco_init
        self.mu = mu
        self.min_step = min_step
        self.window = int(window)
        self._prev_otw: float | None = None
        self._prev_f: float | None = None
        self._prev_count: int | None = None
        self._n = 0
        self._otw_sum = 0.0
        self._coarse_sum = 0.0
        self.history: list[float] = []

    def update(self, otw: float, count: int, f_ref: float, div: int,
               coarse_hz: float = 0.0) -> float:
        if self._prev_count is None:
            # anchor the counter without counting this cycle, so the first
            # window spans exactly ``window`` reference periods
            self._prev_count = count
            return self.kdco
        self._otw_sum += otw
        self._coarse_sum += coarse_hz
        self._n += 1
        if self._n < self.window:
            return self.kdco

        f_now = ((count - self._prev_count) * f_ref * div
                 - self._coarse_sum) / self._n
        otw_now = self._otw_sum / self._n
        self._prev_count, self._n = count, 0
        self._otw_sum = self._coarse_sum = 0.0

        if self._prev_f is not None:
            d_otw = otw_now - self._prev_otw
            if abs(d_otw) >= self.min_step:
                err = (f_now - self._prev_f) - self.kdco * d_otw
                self.kdco += self.mu * err * d_otw / (d_otw ** 2 + 1e-30)
                self.kdco = max(self.kdco, 1.0)
                self.history.append(self.kdco)
        self._prev_f, self._prev_otw = f_now, otw_now
        return self.kdco
