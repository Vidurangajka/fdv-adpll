"""MASH modulators and the shaped-quantisation-noise reference.

The MASH is the architecture the paper argues against: its shaped noise is what
the FDVPD avoids by subtracting the fractional word in the voltage domain
(Sec. II-A).  Getting its noise right is therefore part of the comparison.
"""

import math

import numpy as np
import pytest

from fdvadpll.dsp import welch_psd
from fdvadpll.sdm import MashSdm, mash_sequence, shaped_qn_psd


F_REF = 80e6


@pytest.mark.parametrize("order", [1, 2, 3])
@pytest.mark.parametrize("frac", [0.25, 0.1, 0.2968750])
def test_mean_output_equals_the_input(order, frac):
    """A modulator that does not track its input is useless."""
    y = mash_sequence(frac, 1 << 16, order=order)
    assert y.mean() == pytest.approx(frac, abs=2e-4)


@pytest.mark.parametrize("order", [1, 2, 3])
def test_output_is_integer_valued(order):
    y = mash_sequence(0.3, 4096, order=order)
    assert np.allclose(y, np.round(y))


def test_mash1_output_is_binary():
    y = mash_sequence(0.3, 4096, order=1)
    assert set(np.unique(y)) <= {0.0, 1.0}


@pytest.mark.parametrize("order,expected_range", [(2, 3), (3, 7)])
def test_higher_order_widens_the_output_swing(order, expected_range):
    """MASH-1-1 spans 2 LSB, MASH-1-1-1 spans 4 -- the divider modulus range."""
    y = mash_sequence(0.3, 1 << 14, order=order)
    assert y.max() - y.min() <= expected_range


def test_invalid_order_is_rejected():
    with pytest.raises(ValueError, match="order must be"):
        MashSdm(order=4)


@pytest.mark.parametrize("order,slope", [(2, 40.0), (3, 60.0)])
def test_noise_shaping_slope(order, slope):
    """An Lth-order MASH shapes its quantisation noise as f^(2L).

    Order 1 is left out on purpose: its quantisation error is strongly
    correlated with the input, so the spectrum is a comb of tones rather than a
    shaped noise floor and a slope fit is meaningless (see
    :func:`test_mash1_is_tonal_where_higher_orders_are_not`).
    """
    frac = 0.1234567
    y = mash_sequence(frac, 1 << 17, order=order,
                      rng=np.random.default_rng(0))
    f, s = welch_psd(y - frac, F_REF, nperseg=1 << 13)
    band = (f > 2e4) & (f < 2e6)          # well below f_REF/2, away from DC
    got = 10.0 * np.polyfit(np.log10(f[band]), np.log10(s[band]), 1)[0]
    assert got == pytest.approx(slope, abs=4.0)


def test_mash1_is_tonal_where_higher_orders_are_not():
    """Why an MMDIV-based synthesiser needs at least MASH-1-1."""
    frac = 0.1234567
    peaks = []
    for order in (1, 2):
        y = mash_sequence(frac, 1 << 15, order=order)
        _, s = welch_psd(y - frac, F_REF, nperseg=1 << 12)
        peaks.append(s.max() / np.median(s))
    assert peaks[0] > 10 * peaks[1]


def test_finite_accumulator_width_quantises_the_input():
    """With n_bits set, the mean can only land on a 2^-n grid."""
    m = MashSdm(order=1, n_bits=4)
    y = np.array([m.step(0.3) for _ in range(1 << 15)])
    assert y.mean() == pytest.approx(4.0 / 16.0, abs=2e-3)   # floor(0.3*16)/16


def test_reset_restores_the_initial_state():
    m = MashSdm(order=2)
    first = [m.step(0.3) for _ in range(64)]
    m.reset()
    again = [m.step(0.3) for _ in range(64)]
    assert np.allclose(first, again)


def test_integer_part_passes_through():
    m = MashSdm(order=1)
    y = np.array([m.step(21.25) for _ in range(4096)])
    assert y.mean() == pytest.approx(21.25, abs=1e-3)
    assert set(np.unique(y)) <= {21.0, 22.0}


def test_dither_breaks_up_the_idle_tone():
    """A rational input like 1/4 is the classic MASH-1 tone generator."""
    plain = mash_sequence(0.25, 1 << 14, order=1)
    dithered = MashSdm(order=1, n_bits=8, rng=np.random.default_rng(1))
    y = np.array([dithered.step(0.25) for _ in range(1 << 14)])
    f, s_plain = welch_psd(plain - 0.25, F_REF, nperseg=1 << 11)
    f, s_dith = welch_psd(y - y.mean(), F_REF, nperseg=1 << 11)
    assert s_dith.max() / np.median(s_dith) < s_plain.max() / np.median(s_plain)


# ------------------------------------------------ the analytic comparison --
def test_shaped_qn_psd_follows_the_sine_law():
    f = np.array([1e5, 1e6])
    t_step = 1.0 / 3.36e9
    for order in (1, 2, 3):
        L = shaped_qn_psd(f, F_REF, order, t_step)
        q = (2 * math.pi * t_step * F_REF) ** 2 / 12.0 / F_REF
        want = 0.5 * q * np.abs(2 * np.sin(np.pi * f / F_REF)) ** (2 * order)
        assert np.allclose(L, want)


def test_shaped_qn_rises_with_offset():
    f = np.logspace(4, math.log10(F_REF / 2), 200)
    L = shaped_qn_psd(f, F_REF, 2, 1.0 / 3.36e9)
    assert np.all(np.diff(L) > 0)


def _loop_filtered_jitter(f, L):
    """Jitter a given PD-referred profile contributes at the output."""
    from fdvadpll import default_design, noise

    d = default_design()
    h_lp, _ = noise.loop_transfer(d, f)
    return noise.integrated_jitter(f, L * np.abs(h_lp) ** 2, d.f_ckv)


def test_mmdiv_quantisation_swamps_the_fdvpd_floor():
    """Sec. II-A's argument, in numbers.

    A MASH-dithered multi-modulus divider quantises the phase to one whole CKV
    period, and shapes that error *upwards*; the FDVPD subtracts the fractional
    word in the voltage domain and quantises the residue to ~12 bit.  Deep in
    band the shaped noise is negligible, so the fair comparison is what each
    contributes to the integrated jitter after the loop filter -- where the
    MMDIV alone costs as much jitter as the paper's entire measured 82 fs.
    """
    from fdvadpll import default_design, noise

    d = default_design()
    f = np.logspace(4, math.log10(40e6), 4000)
    mmdiv = _loop_filtered_jitter(f, shaped_qn_psd(f, d.f_ref, 2, d.t_ckv))
    fdvpd = _loop_filtered_jitter(
        f, np.full_like(f, noise.pd_noise_floor(d, 5).quantisation))
    assert fdvpd * 1e15 < 15.0
    assert mmdiv > 4 * fdvpd
    assert mmdiv * 1e15 > 50.0        # same order as the whole measured budget


def test_higher_mash_order_is_worse_behind_a_type2_loop():
    """Steeper shaping loses to a loop that only rolls off at 20 dB/decade.

    Counter-intuitive but real, and the reason an MMDIV synthesiser needs extra
    filter poles that this architecture does not.
    """
    from fdvadpll import default_design

    d = default_design()
    f = np.logspace(4, math.log10(40e6), 4000)
    j = [_loop_filtered_jitter(f, shaped_qn_psd(f, d.f_ref, o, d.t_ckv))
         for o in (1, 2, 3)]
    assert j[0] < j[1] < j[2]
