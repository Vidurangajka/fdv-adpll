"""Digital loop filter and the counter-based frequency-lock path (Fig. 2)."""

from __future__ import annotations

import math

import numpy as np

from .params import DesignParams

__all__ = ["DigitalLoopFilter", "FrequencyLockLoop"]


class DigitalLoopFilter:
    """Normalised type-II proportional-integral loop filter.

        NTW[k] = alpha * phi_e[k] + I[k],      I[k] = I[k-1] + rho * phi_e[k]

    with ``phi_e`` the phase error in *output* (CKV) cycles.  The continuous
    equivalents are ``omega_u = alpha*f_REF``, ``omega_n = sqrt(rho)*f_REF``
    and ``zeta = alpha/(2*sqrt(rho))``, so a target bandwidth and damping map
    directly onto the coefficients (see :meth:`params.LoopParams.alpha_rho`).

    Optional single-pole IIR stages attenuate the reference spur and the
    quantisation noise above the loop bandwidth, as in a conventional ADPLL.
    """

    def __init__(self, design: DesignParams, iir_poles=()):
        self.design = design
        self.loop = design.loop
        self.alpha, self.rho = self.loop.alpha_rho(design.f_ref)
        self.alpha_gear, self.rho_gear = self.loop.alpha_rho(
            design.f_ref, self.loop.gear_bandwidth, self.loop.gear_damping)
        self.integrator = 0.0
        self.k = 0
        self.iir_lambda = np.asarray(iir_poles, dtype=float)
        self.iir_state = np.zeros(len(self.iir_lambda))
        self.freeze_integrator = False

    # ------------------------------------------------------------------
    @property
    def gearing(self) -> bool:
        return self.k < self.loop.gear_cycles

    def coefficients(self):
        if self.gearing:
            return self.alpha_gear, self.rho_gear
        return self.alpha, self.rho

    def bandwidth(self) -> float:
        a, _ = self.coefficients()
        return a * self.design.f_ref / (2.0 * math.pi)

    def reset(self, integrator: float = 0.0) -> None:
        self.integrator = integrator
        self.iir_state[:] = 0.0
        self.k = 0

    def step(self, phi_e: float) -> float:
        """One update.  ``phi_e`` in output cycles, returns the tuning word."""
        alpha, rho = self.coefficients()
        x = phi_e
        for i, lam in enumerate(self.iir_lambda):
            self.iir_state[i] = (1.0 - lam) * self.iir_state[i] + lam * x
            x = self.iir_state[i]
        if not self.freeze_integrator:
            self.integrator += rho * x
        self.k += 1
        return alpha * x + self.integrator


class FrequencyLockLoop:
    """Counter-based frequency/integer-phase path (dashed block of Fig. 2).

    A low-power counter quantises the integer CKV phase; the frequency variance
    is taken from two successive sampled read-outs, compared against FCW to
    form a frequency error, and accumulated into the integer phase error.  It
    removes the "blind" bang-bang period while the residue is outside the 7 b
    ADC range, and is switched off once the ADC has stayed in range for
    ``settle_cycles`` consecutive reference periods (Sec. II-A, Sec. III-A).

    Two details decide whether this works in a fractional channel.

    *The reference side is accumulated and then floored*, rather than compared
    against ``FCW_pd`` directly.  The counter can only ever report a whole
    number of CKVd edges, so a raw ``FCW_pd - d_count`` never reaches zero in a
    fractional channel: it swings by the fractional word every cycle.
    Differencing the floored accumulator gives the integer number of edges that
    genuinely should have elapsed, which does go to zero in lock.

    *The measurement is averaged over ``window`` reference periods.*  One
    reference period resolves frequency only to one CKVd edge, i.e. ``f_REF``
    at the divided clock and ``f_REF*div`` -- 160 MHz -- at the output.  The
    phase detector, by contrast, only captures a few megahertz: its whole range
    is ``adc_full_scale/SR`` out of one PD period.  Correcting on a single
    read-out therefore kicks the oscillator clean out of the detector's range
    every time the quantisation flips, and the loop limit-cycles forever
    instead of settling -- which is exactly what a static lock offset makes it
    do in a fractional channel, while an integer channel sits at a zero error
    and never notices.  Averaging over ``window`` cycles divides the
    quantisation by ``window`` and brings the correction back inside the
    capture range.
    """

    def __init__(self, design: DesignParams, gain_freq: float = 0.25,
                 gain_phase: float = 0.02, settle_cycles: int = 256,
                 window: int = 64):
        self.design = design
        self.fcw_pd = design.fcw_pd
        self.gain_freq = gain_freq
        self.gain_phase = gain_phase
        self.settle_cycles = settle_cycles
        self.window = max(int(window), 1)
        self.prev_count: int | None = None
        self.phase_err = 0.0
        self.enabled = True
        self._in_range = 0
        self.disabled_at: int | None = None
        self._expected = 0.0            # accumulated ideal CKVd phase
        self._prev_expected_int = 0
        self._n = 0
        self._out = 0.0                 # held between measurements

    @property
    def frequency_resolution(self) -> float:
        """Smallest output frequency step this detector can resolve [Hz]."""
        return self.design.f_ref * self.design.fb_div / self.window

    def step(self, count: int, adc_saturated: bool, k: int) -> float:
        """Return the contribution to the tuning word, in output cycles/ref."""
        if not self.enabled:
            return 0.0

        self._in_range = 0 if adc_saturated else self._in_range + 1
        if self._in_range >= self.settle_cycles:
            self.enabled = False
            self.disabled_at = k
            return 0.0

        expected_int = int(self._expected // 1.0)
        self._expected += self.fcw_pd

        if self.prev_count is None:
            self.prev_count = count
            self._prev_expected_int = expected_int
            return 0.0

        self._n += 1
        if self._n < self.window:
            return self._out            # hold the last correction

        d = count - self.prev_count
        d_expected = expected_int - self._prev_expected_int
        self.prev_count = count
        self._prev_expected_int = expected_int
        self._n = 0

        # CKVd cycles per REF cycle, averaged over the window
        freq_err = (d_expected - d) / self.window
        self.phase_err += freq_err
        # referred to output cycles
        div = self.design.fb_div
        self._out = div * (self.gain_freq * freq_err
                           + self.gain_phase * self.phase_err)
        return self._out
