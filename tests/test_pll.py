"""End-to-end behaviour of the simulator.

The load-bearing tests here are the cross-validations: an event-driven run and
the independent s-domain model of :mod:`fdvadpll.noise` must agree.  They share
no code beyond the design parameters, so agreement to a few per cent is real
evidence that both are right.
"""

import math
from dataclasses import replace

import numpy as np
import pytest

from fdvadpll import (FdvPll, default_design, fractional_design,
                      measured_fit_design, noise)


F_GRID = np.logspace(3, math.log10(40e6), 4000)


def _analytic_jitter(design, frac_bit=None):
    L = noise.output_phase_noise(design, F_GRID, frac_bit)
    return noise.integrated_jitter(F_GRID, L, design.f_ckv)


def _settling_cycle(result, threshold=5e-12) -> int:
    """First cycle after which tau stays within ``threshold`` of its final value.

    Measured against the settled value, not against zero: the loop locks with a
    static offset of one V_OS margin (~100 ps), so an absolute test would say
    nothing ever settles.  And not ``argmax(|tau| < threshold)`` either -- tau
    starts at exactly zero, which would report instant settling every time.
    """
    tau = result.tau
    final = float(np.median(tau[-len(tau) // 8:]))
    late = np.nonzero(np.abs(tau - final) > threshold)[0]
    return int(late[-1]) + 1 if len(late) else 0


# ------------------------------------------------------------- basic lock --
def test_integer_channel_locks(integer_run):
    r = integer_run
    assert r.cycle_slips == 0
    assert not r.saturated[r.discard:].any()
    assert np.std(r.tau[r.discard:]) < 1e-12        # sub-picosecond residue


def test_fractional_channel_locks(fractional_run):
    r = fractional_run
    assert r.cycle_slips == 0
    assert not r.saturated[r.discard:].any()


def test_output_frequency_is_the_programmed_one(integer_run):
    r = integer_run
    assert np.mean(r.f_dco[r.discard:]) == pytest.approx(r.f_ckv, rel=1e-5)


def test_fractional_output_frequency(fractional_run):
    r = fractional_run
    assert np.mean(r.f_dco[r.discard:]) == pytest.approx(r.f_ckv, rel=1e-5)
    # ...and it really is a different channel
    assert r.design.fcw != round(r.design.fcw)


def test_runs_are_reproducible():
    a = FdvPll(default_design(), seed=3).run(1 << 12)
    b = FdvPll(default_design(), seed=3).run(1 << 12)
    assert np.array_equal(a.tau, b.tau)


def test_seed_changes_the_realisation():
    a = FdvPll(default_design(), seed=3).run(1 << 12)
    b = FdvPll(default_design(), seed=4).run(1 << 12)
    assert not np.array_equal(a.tau, b.tau)


def test_rejects_an_unknown_integer_path():
    with pytest.raises(ValueError, match="integer_path"):
        FdvPll(default_design(), integer_path="magic")


def test_result_arrays_are_all_one_per_reference_cycle(integer_run):
    n = len(integer_run.tau)
    for name in ("phi_e", "code", "dac_code", "dt", "ntw", "f_dco",
                 "saturated", "n_slip"):
        assert len(getattr(integer_run, name)) == n


def test_meta_records_the_configuration():
    r = FdvPll(default_design(), seed=1, dco_mode="banks",
               integer_path="phase", calibrate=("gain",)).run(1 << 12)
    assert r.meta["dco_mode"] == "banks"
    assert r.meta["integer_path"] == "phase"
    assert r.meta["calibrate"] == ["gain"]
    assert r.meta["seed"] == 1


# ------------------------------------------------- acquisition behaviour --
def test_fll_hands_over_to_the_phase_detector(integer_run):
    """It must switch off, and early -- it is an acquisition aid only."""
    assert integer_run.fll_disabled_at is not None
    assert integer_run.fll_disabled_at < 5000


def test_fll_acquires_a_large_initial_frequency_error():
    r = FdvPll(default_design(), seed=2, f_init_offset=-30e6).run(1 << 15)
    assert r.fll_disabled_at is not None
    assert not r.saturated[r.discard:].any()
    assert np.mean(r.f_dco[r.discard:]) == pytest.approx(r.f_ckv, rel=1e-5)


def test_adc_saturates_before_it_locks():
    """The bang-bang regime the counter path exists to shorten."""
    r = FdvPll(default_design(), seed=2, f_init_offset=-30e6).run(1 << 15)
    assert r.saturated[:200].any()


@pytest.mark.parametrize("path", ["fll", "phase", "none"])
def test_all_integer_paths_lock(path):
    r = FdvPll(default_design(), seed=4, integer_path=path).run(1 << 14)
    assert r.cycle_slips == 0
    assert r.jitter() < 200e-15


def test_phase_path_does_not_use_the_fll():
    r = FdvPll(default_design(), seed=4, integer_path="phase").run(1 << 12)
    assert r.fll_disabled_at is None


def test_gear_shift_speeds_up_acquisition():
    d = default_design()
    slow = replace(d, loop=replace(d.loop, gear_cycles=0))
    a = FdvPll(d, seed=2, f_init_offset=-20e6).run(1 << 14)
    b = FdvPll(slow, seed=2, f_init_offset=-20e6).run(1 << 14)
    assert _settling_cycle(a) < _settling_cycle(b)


# --------------------------------------- cross-validation against Sec. II-B --
@pytest.mark.slow
def test_integer_jitter_matches_the_analytic_model():
    d = default_design()
    r = FdvPll(d, seed=11).run(1 << 17)
    assert r.jitter() == pytest.approx(_analytic_jitter(d), rel=0.10)


@pytest.mark.slow
def test_measured_fit_jitter_matches_the_analytic_model():
    d = measured_fit_design()
    r = FdvPll(d, seed=11).run(1 << 17)
    assert r.jitter() == pytest.approx(_analytic_jitter(d), rel=0.10)
    # ...and therefore lands on the paper's measured 82 fs
    assert r.jitter() * 1e15 == pytest.approx(82.0, abs=8.0)


@pytest.mark.slow
@pytest.mark.parametrize("bit", [4, 6, 8, 10])
def test_fractional_jitter_matches_the_analytic_model(bit):
    """Regression: deep fractional channels used to limit-cycle instead of lock.

    Two separate faults did it -- the FLL comparing an integer count against a
    fractional FCW, and the counter double-counting a cycle whenever the ramp
    window swept through zero.  Both were invisible in an integer channel.
    """
    d = fractional_design(bit, measured_fit_design())
    r = FdvPll(d, seed=5).run(1 << 16)
    assert not r.saturated[r.discard:].any()
    assert r.jitter() == pytest.approx(_analytic_jitter(d, bit), rel=0.12)


@pytest.mark.slow
def test_spectrum_follows_the_analytic_profile():
    """Not just the integral: the shape has to match across the band too."""
    d = measured_fit_design()
    r = FdvPll(d, seed=11).run(1 << 17)
    f, L = r.phase_noise()
    model = noise.output_phase_noise(d, f, None)
    band = (f > 3e4) & (f < 1e7)
    err = 10 * np.log10(L[band] / model[band])
    assert abs(np.median(err)) < 2.0
    assert np.std(err) < 6.0


# ------------------------------------------------------ noise attribution --
@pytest.mark.slow
def test_switching_off_every_noise_source_leaves_almost_nothing():
    r = FdvPll(default_design(), seed=1, dco_noise=False, ref_noise=False,
               pd_noise=False).run(1 << 15)
    # only the DAC/ADC quantisation grid is left
    assert r.jitter() * 1e15 < 15.0


def test_pd_noise_flag_disables_every_thermal_term():
    """Regression: it used to leave kT/C on, i.e. about one ADC LSB."""
    pll = FdvPll(default_design(), pd_noise=False)
    f = pll.design.fdvpd
    assert f.comparator_noise == 0.0
    assert f.gamma == 0.0
    assert not f.ktc_enabled
    assert noise.L_ktc(pll.design) == 0.0


@pytest.mark.slow
@pytest.mark.parametrize("off", ["dco_noise", "ref_noise", "pd_noise"])
def test_each_noise_source_contributes(off):
    full = FdvPll(default_design(), seed=1).run(1 << 15).jitter()
    less = FdvPll(default_design(), seed=1, **{off: False}).run(1 << 15).jitter()
    assert less < full


# ------------------------------------------------------------------ spurs --
def test_integer_channel_reports_no_fractional_spur(integer_run):
    f, dbc = integer_run.fractional_spur()
    assert math.isnan(f)
    assert dbc == float("-inf")


@pytest.mark.slow
def test_ideal_fractional_channel_has_a_low_spur():
    """With a perfect I-DAC the fractional word cancels in the voltage domain."""
    d = fractional_design(5, measured_fit_design())
    r = FdvPll(d, seed=5).run(1 << 16)
    f, dbc = r.fractional_spur()
    assert f == pytest.approx(abs(d.fcw - round(d.fcw)) * d.f_ref, rel=1e-3)
    assert dbc < -70.0


@pytest.mark.slow
def test_dac_gain_error_raises_the_fractional_spur():
    """Sec. III-E's dominant spur mechanism, and what the gain LMS removes."""
    base = fractional_design(5, measured_fit_design())
    bad = replace(base, fdvpd=replace(base.fdvpd, dac_gain_err=0.01))
    a = FdvPll(base, seed=5).run(1 << 16).fractional_spur()[1]
    b = FdvPll(bad, seed=5).run(1 << 16).fractional_spur()[1]
    assert b > a + 15.0


@pytest.mark.slow
def test_gain_calibration_recovers_the_spur():
    base = fractional_design(5, measured_fit_design())
    bad = replace(base, fdvpd=replace(base.fdvpd, dac_gain_err=0.01))
    off = FdvPll(bad, seed=5).run(1 << 16).fractional_spur()[1]
    on = FdvPll(bad, seed=5, calibrate=("gain",),
                cal_mu={"gain": 5e-3}).run(1 << 16)
    assert on.gain_history is not None and len(on.gain_history) > 0
    assert on.fractional_spur()[1] < off


@pytest.mark.slow
def test_dac_mismatch_raises_the_spur():
    base = fractional_design(5, measured_fit_design())
    bad = replace(base, fdvpd=replace(base.fdvpd, dac_sigma_lsb=1.0))
    a = FdvPll(base, seed=5).run(1 << 16)
    b = FdvPll(bad, seed=5).run(1 << 16)
    assert b.meta["dac_inl_lsb"] > a.meta["dac_inl_lsb"]
    assert b.fractional_spur()[1] > a.fractional_spur()[1] + 5.0


def test_worst_spur_returns_a_frequency_in_band(fractional_run):
    f, dbc = fractional_run.worst_spur(f_lo=10e3)
    assert 10e3 <= f <= fractional_run.f_ref / 2
    assert dbc < 0.0


# --------------------------------------------------------- the ADPLL banks --
@pytest.mark.slow
def test_bank_model_locks_and_costs_a_little_jitter():
    """The quantised tuning banks add the dither the ideal DCO does not have."""
    d = measured_fit_design()
    ideal = FdvPll(d, seed=3, dco_mode="ideal").run(1 << 15)
    banks = FdvPll(d, seed=3, dco_mode="banks").run(1 << 15)
    assert banks.cycle_slips == 0
    assert not banks.saturated[banks.discard:].any()
    assert ideal.jitter() < banks.jitter() < 2.0 * ideal.jitter()


@pytest.mark.slow
def test_bank_model_recovers_from_a_coarse_bank_excursion():
    """Regression: acquisition used to strand the coarse bank off-target."""
    r = FdvPll(measured_fit_design(), seed=3, dco_mode="banks",
               f_init_offset=-25e6).run(1 << 15)
    assert not r.saturated[r.discard:].any()
    assert np.mean(r.f_dco[r.discard:]) == pytest.approx(r.f_ckv, rel=1e-5)


# ----------------------------------------------------------- loop filtering --
@pytest.mark.slow
def test_iir_poles_lower_the_far_out_noise():
    d = measured_fit_design()
    plain = FdvPll(d, seed=6).run(1 << 16)
    filt = FdvPll(d, seed=6, iir_poles=(0.3, 0.3)).run(1 << 16)
    f, a = plain.phase_noise()
    _, b = filt.phase_noise()
    band = (f > 5e6) & (f < 3e7)
    assert np.median(b[band]) < np.median(a[band])


@pytest.mark.slow
def test_wider_loop_tracks_a_frequency_step_faster():
    d = measured_fit_design()
    wide = replace(d, loop=replace(d.loop, bandwidth=2e6, gear_cycles=0))
    narrow = replace(d, loop=replace(d.loop, bandwidth=200e3, gear_cycles=0))
    a = FdvPll(wide, seed=2, f_init_offset=-2e6).run(1 << 14)
    b = FdvPll(narrow, seed=2, f_init_offset=-2e6).run(1 << 14)
    assert _settling_cycle(a) < _settling_cycle(b)


# ------------------------------------------------------------- reporting --
def test_summary_reports_the_headline_numbers(fractional_run):
    text = fractional_run.summary()
    for key in ("f_CKV", "rms jitter", "FoM", "worst spur", "cycle slips",
                "fractional spur"):
        assert key in text


def test_fom_uses_the_measured_power(integer_run):
    assert integer_run.power() == pytest.approx(
        integer_run.design.power.p_total_measured)
    assert integer_run.fom() == pytest.approx(
        noise.fom(integer_run.jitter(), integer_run.power()))


def test_jitter_band_can_be_narrowed(integer_run):
    wide = integer_run.jitter(10e3, 40e6)
    narrow = integer_run.jitter(1e6, 40e6)
    assert narrow < wide


def test_tau_settled_skips_the_discard_window(integer_run):
    assert len(integer_run.tau_settled) == len(integer_run.tau) - integer_run.discard
