"""Oscillator model: period jitter, the ideal tuning law, and the ADPLL banks."""

import math
from dataclasses import replace

import numpy as np
import pytest

from fdvadpll import default_design
from fdvadpll.dco import Dco


def _dco(design=None, mode="ideal", noise=False, seed=0, n_edges=1 << 16):
    d = default_design() if design is None else design
    return Dco(d, np.random.default_rng(seed), n_edges, mode=mode, noise=noise)


def test_rejects_an_unknown_mode():
    with pytest.raises(ValueError, match="ideal.*banks"):
        _dco(mode="magic")


# -------------------------------------------------------------- ideal mode --
def test_ideal_tuning_law(design):
    dco = _dco()
    for ntw in (-0.5, 0.0, 0.25):
        dco.set_normalised_tuning(ntw)
        assert dco.f_dco == pytest.approx(design.dco.f_center
                                          + ntw * design.f_ref)


def test_periods_follow_the_tuning_word(design):
    dco = _dco()
    dco.set_normalised_tuning(0.0)
    p = dco.periods(16)
    assert np.allclose(p, design.fb_div / design.dco.f_center)
    assert len(p) == 16


def test_periods_are_referred_to_the_divided_clock(design):
    """The simulator advances CKVd edges, so one period is fb_div CKV periods."""
    dco = _dco()
    dco.set_normalised_tuning(0.0)
    assert dco.periods(1)[0] == pytest.approx(design.fb_div / design.f_ckv)


def test_period_buffer_grows_instead_of_running_out():
    dco = _dco(n_edges=64)
    dco.set_normalised_tuning(0.0)
    total = 0
    for _ in range(50):
        total += len(dco.periods(32))
    assert total == 50 * 32


# -------------------------------------------------------------- the noise --
def test_divider_scales_the_period_jitter(design):
    """Summing fb_div CKV periods multiplies the period variance by fb_div."""
    dco = _dco(noise=True)
    assert dco.sigma_edge == pytest.approx(
        design.dco.period_jitter_white() * math.sqrt(design.fb_div))


def test_injected_jitter_reproduces_the_target_phase_noise(design):
    """A free-running run must measure back the pn_1mhz it was built from."""
    from fdvadpll import dsp

    dco = _dco(noise=True, n_edges=1 << 19)
    dco.set_normalised_tuning(0.0)
    p = dco.periods(1 << 19)
    tau = np.cumsum(p - p.mean())
    f, L = dsp.phase_noise_from_time_error(tau, design.f_pd, design.f_ckv,
                                           nperseg=1 << 15)
    band = (f > 3e6) & (f < 3e7)       # above the flicker corner
    model = dco.free_running_pn(f[band])
    assert np.median(L[band] / model) == pytest.approx(1.0, rel=0.15)


def test_noise_can_be_switched_off():
    dco = _dco(noise=False)
    assert np.all(dco.dt_noise == 0.0)
    dco.set_normalised_tuning(0.0)
    p = dco.periods(1000)
    assert np.ptp(p) == 0.0          # ptp, not std: std leaves a 1e-25 residue


def test_free_running_pn_has_the_right_shape(design):
    dco = _dco()
    f = np.array([1e6, 1e7])
    L = dco.free_running_pn(f)
    assert 10 * math.log10(L[0]) == pytest.approx(design.dco.pn_1mhz, abs=0.5)
    # 1/f^2 far out
    assert 10 * math.log10(L[0] / L[1]) == pytest.approx(20.0, abs=0.5)
    # 1/f^3 below the corner
    lo = dco.free_running_pn(np.array([1e3, 1e4]))
    assert 10 * math.log10(lo[0] / lo[1]) == pytest.approx(30.0, abs=0.5)


# --------------------------------------------------------- the ADPLL banks --
def test_banks_track_the_ideal_law_within_one_lsb(design):
    """The quantised banks must land within a tracking LSB of the ideal curve."""
    dco = _dco(mode="banks")
    for ntw in np.linspace(-0.4, 0.4, 41):
        dco.set_normalised_tuning(float(ntw))
        ideal = design.dco.f_center + ntw * design.f_ref
        assert abs(dco.f_dco - ideal) < 2 * design.dco.kdco


def test_banks_are_monotonic_in_the_tuning_word():
    dco = _dco(mode="banks")
    f = []
    for ntw in np.linspace(-0.4, 0.4, 200):
        dco.set_normalised_tuning(float(ntw))
        f.append(dco.f_dco)
    # allow the sigma-delta dither to wobble by an LSB, but no reversals
    assert np.polyfit(np.arange(len(f)), f, 1)[0] > 0
    assert np.all(np.diff(f) > -2 * 12e3)


def test_coarse_bank_engages_only_on_overflow(design):
    """Inside the tracking range the acquisition bank must stay put."""
    dco = _dco(mode="banks")
    trk_range = 2 ** (design.dco.trk_bits - 1) * design.dco.kdco
    dco.set_normalised_tuning(0.5 * trk_range / design.f_ref)
    assert dco.otw_acq == 0
    dco.set_normalised_tuning(10.0 * trk_range / design.f_ref)
    assert dco.otw_acq != 0


def test_coarse_bank_contribution_is_removed_from_the_tracking_word(design):
    """Regression: a stale coarse offset used to be double-counted.

    Once acquisition moved the coarse bank, every later cycle recomputed the
    total tuning word from scratch; if its coarse share was not subtracted
    again the oscillator sat megahertz off target and the ADC never came out
    of saturation.
    """
    dco = _dco(mode="banks")
    # push hard enough to move the coarse bank...
    dco.set_normalised_tuning(0.4)
    assert dco.otw_acq != 0
    # ...then come back to a word the tracking bank can serve on its own
    dco.set_normalised_tuning(0.0)
    assert dco.f_dco == pytest.approx(design.dco.f_center,
                                      abs=2 * design.dco.kdco)
    assert dco.otw_acq == 0          # the coarse bank released, not stranded
    assert dco.coarse_hz == 0.0


def test_coarse_hz_reports_what_the_coarse_banks_contribute(design):
    dco = _dco(mode="banks")
    dco.set_normalised_tuning(0.3)
    assert dco.coarse_hz == pytest.approx(
        dco.otw_pvt * design.dco.kdco_pvt + dco.otw_acq * design.dco.kdco_acq)
    assert dco.f_dco == pytest.approx(
        design.dco.f_center + dco.coarse_hz
        + round(dco.otw_trk) * design.dco.kdco, abs=design.dco.kdco)


def test_tracking_bank_saturates_gracefully(design):
    """Beyond both banks the frequency clamps instead of wrapping."""
    dco = _dco(mode="banks")
    dco.set_normalised_tuning(50.0)
    hi = dco.f_dco
    dco.set_normalised_tuning(500.0)
    assert dco.f_dco == pytest.approx(hi)
    assert dco.f_dco > design.dco.f_center


def test_dither_gives_sub_lsb_average_resolution(design):
    """The MASH makes the *mean* frequency finer than one K_DCO step."""
    dco = _dco(mode="banks")
    target = 10.5 * design.dco.kdco          # deliberately half an LSB
    ntw = target / design.f_ref
    f = []
    for _ in range(4096):
        dco.set_normalised_tuning(ntw)
        f.append(dco.f_dco)
    assert np.mean(f) == pytest.approx(design.dco.f_center + target,
                                       abs=0.1 * design.dco.kdco)
    assert len(set(f)) > 1                   # it really is dithering


def test_kdco_error_only_moves_the_estimate(design):
    """The true gain is a device property; kdco_err perturbs what the loop assumes."""
    d = replace(design, dco=replace(design.dco, kdco_err=0.3))
    dco = _dco(d, mode="banks")
    assert dco.kdco_true == pytest.approx(design.dco.kdco)
    assert dco.kdco_hat == pytest.approx(design.dco.kdco * 1.3)


def test_kdco_error_scales_the_applied_tuning(design):
    """An over-estimated K_DCO under-drives the oscillator, i.e. lowers loop gain."""
    d = replace(design, dco=replace(design.dco, kdco_err=0.3))
    good, bad = _dco(mode="banks"), _dco(d, mode="banks")
    ntw = 0.01
    good.set_normalised_tuning(ntw)
    bad.set_normalised_tuning(ntw)
    df_good = good.f_dco - design.dco.f_center
    df_bad = bad.f_dco - design.dco.f_center
    assert df_bad == pytest.approx(df_good / 1.3, rel=0.05)
