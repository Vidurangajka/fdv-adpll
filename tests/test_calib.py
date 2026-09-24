"""Background calibration loops (the ADPLL extension).

The gain and INL loops are exercised against the mechanism they exist to
cancel; the K_DCO loop is checked open-loop, where its excitation is known,
because in lock the tuning word barely moves and there is nothing to estimate
from (see the class docstring).
"""

import numpy as np
import pytest

from fdvadpll.calib import DacInlCalibration, KdcoCalibration, PdGainCalibration


# ------------------------------------------------------------- PD gain LMS --
def test_gain_correction_starts_transparent():
    cal = PdGainCalibration()
    for x in (0.0, 0.25, 0.9):
        assert cal.correct(x) == pytest.approx(x)


def test_gain_correction_stays_inside_one_period():
    cal = PdGainCalibration(g_init=1.2)
    assert 0.0 <= cal.correct(0.99) < 1.0
    assert cal.correct(-0.5) == 0.0


def test_gain_lms_converges_on_a_known_mismatch():
    """Plant: the measured phase error is (g_true*g_hat - 1) * (T_frac - 1/2)."""
    g_true = 1.0 / 1.01                      # 1 % DAC-to-ramp gain error
    cal = PdGainCalibration(mu=0.05, leak=0.0)
    rng = np.random.default_rng(0)
    for _ in range(20000):
        t_frac = float(rng.uniform(0.0, 1.0))
        phi_e = (1.0 - g_true * cal.g) * (t_frac - 0.5)
        cal.update(phi_e, t_frac)
    assert cal.g == pytest.approx(1.01, rel=0.02)


def test_gain_clamp_cannot_exceed_what_the_detector_absorbs():
    """A wider clamp would let a bad estimate rail the ADC and stop recovery."""
    from fdvadpll import default_design

    d = default_design()
    absorbable = d.fdvpd.adc_full_scale / d.fdvpd.slew_rate / d.t_pd
    assert PdGainCalibration().clamp <= absorbable


def test_leak_pulls_a_starved_estimate_back_to_unity():
    """The recovery path: relax() runs on railed cycles, where update() cannot."""
    cal = PdGainCalibration(leak=1e-3)
    cal.g = 1.03
    for _ in range(20000):
        cal.relax()
    assert cal.g == pytest.approx(1.0, abs=1e-3)


def test_leak_does_not_meaningfully_bias_a_real_gradient():
    g_true = 1.0 / 1.01
    leaky = PdGainCalibration(mu=0.05)
    clean = PdGainCalibration(mu=0.05, leak=0.0)
    rng = np.random.default_rng(0)
    for _ in range(20000):
        t_frac = float(rng.uniform(0.0, 1.0))
        for cal in (leaky, clean):
            cal.update((1.0 - g_true * cal.g) * (t_frac - 0.5), t_frac)
    assert leaky.g == pytest.approx(clean.g, abs=2e-3)


def test_gain_lms_is_clamped():
    cal = PdGainCalibration(mu=0.5, clamp=0.1)
    for _ in range(1000):
        cal.update(10.0, 1.0)
    assert cal.g <= 1.1 + 1e-12
    for _ in range(4000):
        cal.update(-10.0, 1.0)
    assert cal.g >= 0.9 - 1e-12


def test_gain_lms_ignores_an_error_uncorrelated_with_the_fractional_word():
    """White measurement noise must not bias the estimate."""
    cal = PdGainCalibration(mu=1e-3)
    rng = np.random.default_rng(1)
    for _ in range(50000):
        cal.update(float(rng.standard_normal()) * 1e-3,
                   float(rng.uniform(0.0, 1.0)))
    assert cal.g == pytest.approx(1.0, abs=0.02)


def test_gain_history_is_recorded():
    cal = PdGainCalibration()
    for _ in range(10):
        cal.update(0.1, 0.8)
    assert len(cal.history) == 10
    assert cal.history[-1] == cal.g


def test_sign_sign_variant_moves_in_the_same_direction():
    a = PdGainCalibration(mu=1e-3)
    b = PdGainCalibration(mu=1e-3, sign_sign=True)
    for _ in range(100):
        a.update(0.2, 0.9)
        b.update(0.2, 0.9)
    assert (a.g - 1.0) * (b.g - 1.0) > 0


# ------------------------------------------------------------- INL look-up --
def test_inl_lut_indexes_by_thermometer_segment():
    cal = DacInlCalibration(n_segments=16, dac_bits=10, bin_bits=6)
    assert cal.segment(0) == 0
    assert cal.segment(63) == 0
    assert cal.segment(64) == 1
    assert cal.segment(1023) == 15


def test_inl_correction_is_in_dac_lsbs():
    cal = DacInlCalibration(n_segments=16, dac_bits=10, bin_bits=6)
    cal.lut[2] = 4.0
    corrected = cal.correct(0.5, dac_code=2 * 64)
    assert corrected == pytest.approx(0.5 + 4.0 / 1024)


def test_inl_correction_is_clipped_to_one_period():
    cal = DacInlCalibration(n_segments=16, dac_bits=10, bin_bits=6)
    cal.lut[0] = 1e6
    assert cal.correct(0.5, 0) < 1.0


def test_inl_lut_stays_zero_mean():
    """A common-mode entry would be a static phase offset, which the loop owns."""
    cal = DacInlCalibration(n_segments=16, mu=1e-3)
    rng = np.random.default_rng(0)
    for _ in range(5000):
        cal.update(float(rng.standard_normal()), int(rng.integers(0, 1024)))
    assert cal.lut.mean() == pytest.approx(0.0, abs=1e-12)


def test_inl_lms_learns_a_per_segment_offset():
    """Plant: each segment adds a fixed error the LUT should learn to cancel."""
    n_seg = 16
    truth = np.linspace(-3.0, 3.0, n_seg)
    truth -= truth.mean()
    cal = DacInlCalibration(n_segments=n_seg, mu=2e-3, dac_bits=10, bin_bits=6)
    rng = np.random.default_rng(2)
    for _ in range(200000):
        code = int(rng.integers(0, 1024))
        s = code >> 6
        # residual error after the current correction, in DAC LSB.  The sign is
        # the loop's: a positive LUT entry pushes the measured error negative,
        # which is what makes the sign-LMS settle rather than run away.
        phi_e = -(truth[s] + cal.lut[s])
        cal.update(phi_e, code)
    assert np.corrcoef(cal.lut, -truth)[0, 1] > 0.95


# ------------------------------------------------------------ K_DCO tracker --
def _drive_kdco(cal, kdco_true, n_windows, f_ref=80e6, div=2,
                f_center=3.36e9, coarse=0.0, seed=0):
    """Open-loop excitation: a fresh random tuning word every window."""
    rng = np.random.default_rng(seed)
    phase, otw = 0.0, 0.0
    for k in range(n_windows * cal.window + 1):
        if k % cal.window == 0:
            otw = float(rng.integers(-110, 110))
        phase += (f_center + coarse + kdco_true * otw) / f_ref / div
        cal.update(otw, int(phase), f_ref, div, coarse)
    return cal.kdco


def test_kdco_converges_from_a_thirty_percent_error():
    cal = KdcoCalibration(kdco_init=15600.0, mu=0.25, window=16384)
    assert _drive_kdco(cal, 12000.0, 30) == pytest.approx(12000.0, rel=0.02)
    assert len(cal.history) > 10


def test_kdco_converges_from_an_underestimate():
    cal = KdcoCalibration(kdco_init=6000.0, mu=0.25, window=16384)
    assert _drive_kdco(cal, 12000.0, 30) == pytest.approx(12000.0, rel=0.02)


def test_kdco_subtracts_the_coarse_banks():
    """Regression: a coarse-bank offset used to be credited to K_DCO."""
    cal = KdcoCalibration(kdco_init=15600.0, mu=0.25, window=16384)
    got = _drive_kdco(cal, 12000.0, 30, coarse=31.5e6)
    assert got == pytest.approx(12000.0, rel=0.02)


def test_kdco_needs_a_window_long_enough_to_see_a_bank_step():
    """A short window quantises the frequency reading coarser than the bank.

    Every measurement then returns the same integer count, the apparent step is
    zero, and a naive LMS would conclude K_DCO is zero.  The guard is the
    window length, so this documents where the floor is.
    """
    short = KdcoCalibration(kdco_init=12000.0, mu=0.25, window=128)
    got = _drive_kdco(short, 12000.0, 30)
    # resolution is 1.25 MHz against a 3.1 MHz bank: the estimate is not merely
    # noisy, it is meaningless, and it starts from the right answer
    assert abs(got - 12000.0) / 12000.0 > 0.5

    long = KdcoCalibration(kdco_init=12000.0, mu=0.25, window=16384)
    assert _drive_kdco(long, 12000.0, 30) == pytest.approx(12000.0, rel=0.02)


def test_kdco_ignores_steps_below_min_step():
    """No excitation, no update -- the in-lock case."""
    cal = KdcoCalibration(kdco_init=15600.0, mu=0.25, window=1024,
                          min_step=0.5)
    phase = 0.0
    for k in range(20 * 1024 + 1):
        phase += (3.36e9 + 12000.0 * 3.0) / 80e6 / 2     # tuning word frozen
        cal.update(3.0, int(phase), 80e6, 2, 0.0)
    assert cal.history == []
    assert cal.kdco == 15600.0


def test_kdco_never_goes_non_positive():
    cal = KdcoCalibration(kdco_init=12000.0, mu=5.0, window=256)
    _drive_kdco(cal, 12000.0, 40)
    assert cal.kdco >= 1.0
