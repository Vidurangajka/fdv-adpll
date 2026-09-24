"""Noise synthesis and spectral estimation.

These are the foundations everything else is measured with, so they are checked
against closed forms rather than against each other wherever possible.
"""

import math

import numpy as np
import pytest

from fdvadpll import dsp


FS = 80e6


# ------------------------------------------------------------- synthesis --
def test_colored_noise_reproduces_a_white_psd(rng):
    """A flat target PSD must come back flat, at the right level."""
    s0 = 1e-24
    x = dsp.colored_noise(1 << 16, FS, lambda f: np.full_like(f, s0), rng)
    f, s = dsp.welch_psd(x, FS, nperseg=1 << 12)
    assert np.median(s) == pytest.approx(s0, rel=0.05)
    # variance = integral S df over [0, fs/2]
    assert np.var(x) == pytest.approx(s0 * FS / 2.0, rel=0.05)


def test_colored_noise_reproduces_a_slope(rng):
    """1/f^2 in, 1/f^2 out: -20 dB per decade."""
    x = dsp.colored_noise(1 << 17, FS, lambda f: 1e-18 / f ** 2, rng)
    f, s = dsp.welch_psd(x, FS, nperseg=1 << 13)
    band = (f > 3e3) & (f < 3e6)
    slope = np.polyfit(np.log10(f[band]), np.log10(s[band]), 1)[0]
    assert slope == pytest.approx(-2.0, abs=0.12)


def test_colored_noise_flattens_below_f_min(rng):
    """Divergent laws are clamped, so the sequence stays finite."""
    x = dsp.colored_noise(1 << 14, FS, lambda f: 1e-18 / f ** 3, rng,
                          f_min=1e5)
    assert np.all(np.isfinite(x))


def test_one_over_f_slope(rng):
    x = dsp.one_over_f(1 << 17, FS, 1e-18, rng)
    f, s = dsp.welch_psd(x, FS, nperseg=1 << 13)
    band = (f > 1e4) & (f < 1e7)
    slope = np.polyfit(np.log10(f[band]), np.log10(s[band]), 1)[0]
    assert slope == pytest.approx(-1.0, abs=0.12)


def test_time_error_from_phase_noise_round_trip(rng):
    """Synthesise a known L(f), measure it back."""
    target = 1e-13                      # linear, 1/Hz
    f0 = 3.36e9
    x = dsp.time_error_from_phase_noise(
        1 << 16, FS, f0, lambda f: np.full_like(f, target), rng)
    f, L = dsp.phase_noise_from_time_error(x, FS, f0, nperseg=1 << 12)
    assert np.median(L) == pytest.approx(target, rel=0.05)


def test_dco_period_jitter_gives_the_expected_one_over_f_squared(rng):
    """White period jitter sigma -> L(f) = f0^3 sigma^2 / f^2."""
    f0, sigma = 1.68e9, 20e-15
    dt = dsp.dco_period_jitter(1 << 18, f0, sigma, 0.0, rng)
    assert np.std(dt) == pytest.approx(sigma, rel=0.02)
    # accumulate the period deviations into an edge-timing error
    tau = np.cumsum(dt)
    f, L = dsp.phase_noise_from_time_error(tau, f0, f0, nperseg=1 << 14)
    band = (f > 1e6) & (f < 1e8)
    predicted = f0 ** 3 * sigma ** 2 / f[band] ** 2
    ratio = np.median(L[band] / predicted)
    assert ratio == pytest.approx(1.0, rel=0.1)


def test_dco_period_jitter_flicker_crosses_at_the_corner(rng):
    """The 1/f^3 part must equal the 1/f^2 part exactly at f_corner."""
    f0, sigma, fc = 1.68e9, 20e-15, 1e6
    dt = dsp.dco_period_jitter(1 << 19, f0, sigma, fc, rng)
    tau = np.cumsum(dt)
    f, L = dsp.phase_noise_from_time_error(tau, f0, f0, nperseg=1 << 15)
    white = f0 ** 3 * sigma ** 2 / f ** 2
    near = np.abs(f - fc) < 0.15 * fc
    # total = white * (1 + fc/f) -> exactly 2x the white part at f = fc
    assert np.median(L[near] / white[near]) == pytest.approx(2.0, rel=0.15)


# ----------------------------------------------------------- measurement --
def test_welch_psd_drops_dc(rng):
    x = rng.standard_normal(4096) + 5.0
    f, s = dsp.welch_psd(x, FS, nperseg=1024)
    assert f[0] > 0.0
    assert np.median(s) == pytest.approx(2.0 / FS, rel=0.15)


def test_rms_jitter_matches_the_closed_form():
    """L = A/f^2 integrates to 2A(1/f_lo - 1/f_hi) rad^2."""
    f = np.logspace(3, 8, 20001)
    A, f0 = 1e-2, 3.36e9
    L = A / f ** 2
    f_lo, f_hi = 10e3, 40e6
    expected = math.sqrt(2.0 * A * (1.0 / f_lo - 1.0 / f_hi)) / (2 * math.pi * f0)
    assert dsp.rms_jitter(f, L, f0, f_lo, f_hi) == pytest.approx(expected, rel=1e-3)


def test_integrate_phase_noise_of_a_flat_profile():
    f = np.linspace(1e3, 1e6, 100001)
    L = np.full_like(f, 1e-14)
    ipn = dsp.integrate_phase_noise(f, L, 1e3, 1e6)
    assert ipn == pytest.approx(2.0 * 1e-14 * (1e6 - 1e3), rel=1e-6)


def test_integrate_phase_noise_rejects_an_empty_band():
    f = np.array([1e3, 1e4, 1e5])
    with pytest.raises(ValueError, match="spectral points"):
        dsp.integrate_phase_noise(f, np.ones(3), 2e5, 3e5)


# ------------------------------------------------------------------ spurs --
def test_find_spurs_reads_a_known_tone_at_its_true_level():
    """A phase tone of peak deviation phi_pk reads (phi_pk/2)^2 in dBc."""
    n, f0 = 1 << 14, 3.36e9
    f_spur = 200 * FS / n           # exactly on a bin: no scalloping loss
    phi_pk = 1e-3                   # rad
    k = np.arange(n)
    tau = phi_pk / (2 * math.pi * f0) * np.sin(2 * math.pi * f_spur * k / FS)
    f, dbc = dsp.find_spurs(tau, FS, f0)
    expected = 20.0 * math.log10(phi_pk / 2.0)
    i = int(np.argmin(np.abs(f - f_spur)))
    assert dbc[i] == pytest.approx(expected, abs=0.1)


def test_find_spurs_blanks_the_dc_bins():
    tau = np.full(1024, 1e-12)
    f, dbc = dsp.find_spurs(tau, FS, 3.36e9, exclude_dc_bins=4)
    assert np.all(np.isneginf(dbc[:3]))


def test_spur_level_locks_onto_the_nearest_tone():
    n, f0, f_spur = 1 << 14, 3.36e9, 2.5e6
    k = np.arange(n)
    tau = 3e-16 * np.sin(2 * math.pi * f_spur * k / FS)
    # ask slightly off-frequency; the search window must still find it
    f_found, dbc = dsp.spur_level(tau, FS, f0, f_spur * 1.0002)
    assert f_found == pytest.approx(f_spur, rel=1e-3)
    assert dbc == pytest.approx(20 * math.log10(2 * math.pi * f0 * 3e-16 / 2),
                                abs=0.2)
