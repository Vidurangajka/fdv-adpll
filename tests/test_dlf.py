"""Digital loop filter and the counter-based frequency-lock path."""

import math
from dataclasses import replace

import numpy as np
import pytest

from fdvadpll import DesignParams, default_design
from fdvadpll.dlf import DigitalLoopFilter, FrequencyLockLoop


# ------------------------------------------------------------ loop filter --
def test_proportional_and_integral_paths(design):
    """NTW = alpha*phi_e + sum(rho*phi_e), with gearing bypassed."""
    d = replace(design, loop=replace(design.loop, gear_cycles=0))
    dlf = DigitalLoopFilter(d)
    alpha, rho = dlf.alpha, dlf.rho
    assert dlf.step(1.0) == pytest.approx(alpha + rho)
    assert dlf.step(1.0) == pytest.approx(alpha + 2 * rho)
    assert dlf.step(0.0) == pytest.approx(2 * rho)      # integrator holds


def test_bandwidth_matches_the_design_target(design):
    d = replace(design, loop=replace(design.loop, gear_cycles=0))
    assert DigitalLoopFilter(d).bandwidth() == pytest.approx(
        design.loop.bandwidth, rel=1e-9)


def test_gear_shift_starts_wide_then_narrows(design):
    dlf = DigitalLoopFilter(design)
    assert dlf.gearing
    assert dlf.bandwidth() == pytest.approx(design.loop.gear_bandwidth)
    for _ in range(design.loop.gear_cycles):
        dlf.step(0.0)
    assert not dlf.gearing
    assert dlf.bandwidth() == pytest.approx(design.loop.bandwidth)


def test_integrator_freeze(design):
    d = replace(design, loop=replace(design.loop, gear_cycles=0))
    dlf = DigitalLoopFilter(d)
    dlf.step(1.0)
    held = dlf.integrator
    dlf.freeze_integrator = True
    dlf.step(1.0)
    assert dlf.integrator == pytest.approx(held)


def test_reset_clears_the_state(design):
    dlf = DigitalLoopFilter(design)
    for _ in range(100):
        dlf.step(1.0)
    dlf.reset()
    assert dlf.integrator == 0.0
    assert dlf.k == 0
    assert np.all(dlf.iir_state == 0.0)


def test_reset_can_preload_the_integrator(design):
    dlf = DigitalLoopFilter(design)
    dlf.reset(integrator=0.25)
    assert dlf.integrator == 0.25


def test_iir_poles_attenuate_high_frequency_error(design):
    """Extra poles trade settling for out-of-band rejection."""
    d = replace(design, loop=replace(design.loop, gear_cycles=0))
    plain = DigitalLoopFilter(d)
    filtered = DigitalLoopFilter(d, iir_poles=(0.1, 0.1))
    alt = np.where(np.arange(256) % 2 == 0, 1.0, -1.0)   # Nyquist-rate error
    a = np.array([plain.step(x) for x in alt])
    b = np.array([filtered.step(x) for x in alt])
    assert np.std(b[64:]) < 0.2 * np.std(a[64:])


def test_iir_poles_settle_to_the_same_dc_gain(design):
    d = replace(design, loop=replace(design.loop, gear_cycles=0))
    plain = DigitalLoopFilter(d)
    filtered = DigitalLoopFilter(d, iir_poles=(0.2,))
    for _ in range(4000):
        a = plain.step(1.0)
        b = filtered.step(1.0)
    assert b == pytest.approx(a, rel=1e-3)


@pytest.mark.parametrize("zeta", [0.5, 1.0, 2.0])
def test_damping_sets_the_alpha_rho_ratio(design, zeta):
    d = replace(design, loop=replace(design.loop, damping=zeta, gear_cycles=0))
    dlf = DigitalLoopFilter(d)
    assert dlf.alpha / (2 * math.sqrt(dlf.rho)) == pytest.approx(zeta)


def test_type2_loop_nulls_a_frequency_offset(design):
    """The integrator is what makes the static phase error zero."""
    d = replace(design, loop=replace(design.loop, gear_cycles=0))
    dlf = DigitalLoopFilter(d)
    # emulate a first-order plant: phase error accumulates the tuning shortfall
    target, phi, ntw = 1e-3, 0.0, 0.0
    for _ in range(20000):
        phi += target - ntw
        ntw = dlf.step(phi)
    assert ntw == pytest.approx(target, rel=1e-3)
    assert abs(phi) < 1e-3


# ---------------------------------------------------------------- the FLL --
def test_fll_is_silent_on_its_first_call(design):
    """It needs two counter read-outs before it can measure a frequency."""
    fll = FrequencyLockLoop(design)
    assert fll.step(0, True, 0) == 0.0


def test_fll_signs_a_slow_oscillator_positive(design):
    """Counting fewer CKVd edges than FCW means the DCO must speed up."""
    fll = FrequencyLockLoop(design, window=1)
    fll.step(0, True, 0)
    out = fll.step(20, True, 1)          # 20 edges vs FCW_pd = 21
    assert out > 0


def test_fll_signs_a_fast_oscillator_negative(design):
    fll = FrequencyLockLoop(design, window=1)
    fll.step(0, True, 0)
    assert fll.step(22, True, 1) < 0


def test_fll_output_is_referred_to_output_cycles(design):
    """The div-2 must be undone before the word reaches the loop."""
    fll = FrequencyLockLoop(design, gain_freq=1.0, gain_phase=0.0, window=1)
    fll.step(0, True, 0)
    out = fll.step(20, True, 1)
    assert out == pytest.approx(design.fb_div * 1.0)


def test_fll_accumulates_a_phase_error(design):
    fll = FrequencyLockLoop(design, gain_freq=0.0, gain_phase=1.0, window=1)
    fll.step(0, True, 0)
    a = fll.step(20, True, 1)
    b = fll.step(40, True, 2)
    assert b == pytest.approx(2 * a)     # two cycles of a constant 1-edge error


# --- the averaging that makes fractional channels work --------------------
def test_fll_holds_its_output_between_measurements(design):
    """One correction per window, held in between -- not a fresh kick each cycle."""
    fll = FrequencyLockLoop(design, window=8)
    fll.step(0, True, 0)
    outs = [fll.step(20 * k, True, k) for k in range(1, 25)]
    assert len(set(outs)) <= 4           # ~3 windows over 24 cycles
    assert outs[0] == outs[1] == outs[2]


def test_fll_averages_the_counter_quantisation(design):
    """A half-count average must come out half, not alternate between 0 and 1."""
    fll = FrequencyLockLoop(design, gain_freq=1.0, gain_phase=0.0, window=8)
    fll.step(0, True, 0)
    # deliver 8 edges short over 8 cycles -> 1 edge/cycle average
    count, out = 0, 0.0
    for k in range(1, 9):
        count += 20
        out = fll.step(count, True, k)
    assert out == pytest.approx(design.fb_div * 1.0)


def test_fll_frequency_resolution_scales_with_the_window(design):
    """The quantisation floor the window exists to beat."""
    assert FrequencyLockLoop(design, window=1).frequency_resolution == \
        pytest.approx(design.f_ref * design.fb_div)
    assert FrequencyLockLoop(design, window=64).frequency_resolution == \
        pytest.approx(design.f_ref * design.fb_div / 64)


def test_fll_default_window_resolves_finer_than_the_pd_capture_range(design):
    """Otherwise a single correction throws the loop out of the detector.

    The phase detector captures ``adc_full_scale/SR`` of timing error per
    reference period; expressed as a frequency that is the largest step the
    FLL may take without blinding the detector it is trying to hand over to.
    """
    fll = FrequencyLockLoop(design)
    capture_time = design.fdvpd.adc_full_scale / design.fdvpd.slew_rate
    capture_hz = capture_time * design.f_ref * design.f_ckv
    assert fll.frequency_resolution < capture_hz


def test_fll_error_goes_to_zero_in_a_fractional_channel(design):
    """Regression: comparing an integer count against a fractional FCW.

    ``FCW_pd - d_count`` can never be zero when FCW is fractional, so the old
    formulation dithered the tuning word by a full count every cycle and the
    loop never settled.  Accumulating and flooring the reference side fixes it.
    """
    from dataclasses import replace

    d = replace(design, fcw=design.fcw - 2.0 ** -5 * design.fb_div)
    fll = FrequencyLockLoop(d, window=1)
    phase = 0.0
    outs = []
    for k in range(400):
        # a perfectly locked oscillator; the accumulator advances at the end
        # of the cycle, as it does in FdvPll.run
        outs.append(fll.step(int(phase), True, k))
        phase += d.fcw_pd
    assert max(abs(o) for o in outs[50:]) < 1e-9


def test_fll_switches_off_once_the_adc_stays_in_range(design):
    fll = FrequencyLockLoop(design, settle_cycles=16)
    count = 0
    for k in range(64):
        count += 21
        out = fll.step(count, False, k)
        if not fll.enabled:
            break
    assert not fll.enabled
    assert fll.disabled_at == 15        # the 16th consecutive in-range sample
    assert out == 0.0
    assert fll.step(count, False, 99) == 0.0     # and stays off


def test_saturation_restarts_the_settle_counter(design):
    fll = FrequencyLockLoop(design, settle_cycles=8)
    count = 0
    for k in range(7):
        count += 21
        fll.step(count, False, k)
    count += 21
    fll.step(count, True, 7)                     # one saturated sample
    assert fll.enabled
    for k in range(8, 20):
        count += 21
        fll.step(count, False, k)
        if not fll.enabled:
            break
    assert fll.disabled_at == 15                 # 8 clean cycles after the hit


def test_fll_drives_a_frequency_error_to_zero(design):
    """Closed loop on a toy integrator plant."""
    fll = FrequencyLockLoop(design, settle_cycles=10 ** 9)
    # a real counter reads floor(accumulated phase), so the count difference
    # averages the true edge rate instead of rounding it away
    edges_per_ref, phase = 18.0, 0.0
    for k in range(3000):
        phase += edges_per_ref
        edges_per_ref += 0.05 * fll.step(int(phase), True, k) / design.fb_div
    assert edges_per_ref == pytest.approx(design.fcw_pd, abs=0.02)
