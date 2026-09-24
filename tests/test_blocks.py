"""I-DAC, ramp generator, SAR ADC and the assembled phase detector.

The imperfections listed in Sec. III-E are what turn a fractional word into a
spur, so each one is checked in isolation here before :mod:`test_pll` looks for
its effect at the output.
"""

import math
from dataclasses import replace

import numpy as np
import pytest

from fdvadpll import DesignParams, default_design
from fdvadpll.blocks import CurrentDAC, Fdvpd, RampGenerator, SarAdc
from fdvadpll.params import KT, FdvpdParams


def _rng(seed=0):
    return np.random.default_rng(seed)


# ------------------------------------------------------------------ I-DAC --
def test_ideal_dac_is_perfectly_linear():
    dac = CurrentDAC(FdvpdParams(), _rng())
    assert dac.n_codes == 1024
    assert np.max(np.abs(dac.inl)) < 1e-9
    assert np.max(np.abs(dac.dnl)) < 1e-9
    assert np.all(np.diff(dac.levels) > 0)


def test_dac_maps_codes_onto_the_full_pd_period():
    dac = CurrentDAC(FdvpdParams(), _rng())
    assert dac.frac[0] == pytest.approx(0.0)
    assert dac.frac[-1] == pytest.approx(1.0 - 1.0 / dac.n_codes)
    assert np.all(dac.frac < 1.0)


def test_dac_segmentation_weights():
    """4 b thermometer of weight 64 LSB each, plus a 6 b binary array."""
    dac = CurrentDAC(FdvpdParams(), _rng())
    assert len(dac.w_thermo) == 15
    assert np.allclose(dac.w_thermo, 64.0)
    assert np.allclose(dac.w_bin, [1, 2, 4, 8, 16, 32])


def test_dac_gain_error_is_a_pure_scaling():
    a = CurrentDAC(FdvpdParams(), _rng())
    b = CurrentDAC(FdvpdParams(dac_gain_err=0.01), _rng())
    assert np.allclose(b.levels, a.levels * 1.01)
    assert np.max(np.abs(b.inl)) < 1e-9      # gain alone is not non-linearity


def test_unit_element_mismatch_follows_the_sqrt_w_law():
    """A weight-64 element must carry 8x the sigma of the 1 LSB element."""
    sig = 0.05
    w_thermo, w_bin = [], []
    for seed in range(400):
        dac = CurrentDAC(FdvpdParams(dac_sigma_lsb=sig), _rng(seed))
        w_thermo.append(dac.w_thermo)
        w_bin.append(dac.w_bin[0])
    s_thermo = np.std(np.concatenate(w_thermo))
    s_unit = np.std(w_bin)
    assert s_thermo == pytest.approx(sig * 8.0, rel=0.15)
    assert s_unit == pytest.approx(sig * 1.0, rel=0.15)


def test_mismatch_produces_inl_peaking_at_segment_boundaries():
    """The signature the segmented topology is chosen to keep small."""
    dac = CurrentDAC(FdvpdParams(dac_sigma_lsb=0.5), _rng(3))
    assert np.max(np.abs(dac.inl)) > 0.5
    # endpoint fit -> INL pinned to zero at both ends
    assert dac.inl[0] == pytest.approx(0.0, abs=1e-9)
    assert dac.inl[-1] == pytest.approx(0.0, abs=1e-9)


def test_segmentation_beats_a_plain_binary_dac_on_dnl():
    """Sec. III-C's reason for segmenting: the mid-scale DNL jump."""
    p_seg = FdvpdParams(dac_sigma_lsb=0.3, dac_msb_thermo_bits=4)
    p_bin = FdvpdParams(dac_sigma_lsb=0.3, dac_msb_thermo_bits=0)
    dnl_seg = np.mean([np.max(np.abs(CurrentDAC(p_seg, _rng(s)).dnl))
                       for s in range(30)])
    dnl_bin = np.mean([np.max(np.abs(CurrentDAC(p_bin, _rng(s)).dnl))
                       for s in range(30)])
    assert dnl_seg < dnl_bin


# ---------------------------------------------------------------- ramp --
def test_slew_rate_is_two_i_r_over_c_sar():
    p = FdvpdParams()
    ramp = RampGenerator(p)
    assert ramp.sr == pytest.approx(2.0 * p.i_ramp / p.c_sar)
    assert ramp.excursion(100e-12) == pytest.approx(80e-3)   # 0.8 mV/ps


def test_ramp_is_linear_without_distortion_terms():
    ramp = RampGenerator(FdvpdParams())
    dt = np.linspace(0, 500e-12, 51)
    v = ramp.excursion(dt)
    assert np.allclose(v, ramp.sr * dt)


def test_ramp_nonlinearity_bends_the_characteristic():
    ramp = RampGenerator(FdvpdParams(ramp_nl2=0.5))
    dt = 300e-12
    lin = ramp.sr * dt
    assert ramp.excursion(dt) == pytest.approx(lin * (1 + 0.5 * lin))


def test_ramp_mismatch_is_a_gain_error_not_an_asymmetry():
    """Sec. III-A: both I_R sources carry the same lag, so it cannot skew."""
    ramp = RampGenerator(FdvpdParams(ramp_mismatch=0.02))
    assert ramp.sr == pytest.approx(0.8e-3 / 1e-12 * 1.02)
    dt = np.array([-200e-12, 200e-12])
    v = ramp.excursion(dt)
    assert v[0] == pytest.approx(-v[1])       # odd symmetry preserved


def test_ktc_sigma_is_two_kt_over_c():
    """Equation (9): the differential pair samples 2kT/C."""
    p = FdvpdParams()
    assert RampGenerator(p).ktc_sigma == pytest.approx(
        math.sqrt(2.0 * KT / p.c_sar))


def test_ramp_current_noise_reproduces_equation_7():
    """The behavioural sigma must agree with the analytic L_ramp_current."""
    from fdvadpll import noise
    d = default_design()
    ramp = RampGenerator(d.fdvpd)
    t_on = 100e-12

    # sampled voltage -> time -> output phase, once per reference cycle
    sigma_v = ramp.current_noise_sigma(t_on)
    sigma_t = sigma_v / d.fdvpd.slew_rate
    sigma_phi = 2 * math.pi * d.f_ckv * sigma_t
    L = sigma_phi ** 2 / d.f_ref                    # white, one-sided, /2 for SSB
    analytic = noise.L_ramp_current(d, None, include_div2=True)
    assert L / 2.0 == pytest.approx(analytic, rel=1e-9)


def test_ramp_current_noise_grows_as_sqrt_of_the_window():
    ramp = RampGenerator(FdvpdParams())
    a = ramp.current_noise_sigma(100e-12)
    b = ramp.current_noise_sigma(400e-12)
    assert b / a == pytest.approx(2.0)
    assert ramp.current_noise_sigma(0.0) == 0.0


# ------------------------------------------------------------- SAR ADC --
def test_sar_transfer_is_a_midtread_quantiser():
    adc = SarAdc(FdvpdParams(comparator_noise=0.0), _rng())
    assert adc.n_codes == 128
    assert adc.full_scale == pytest.approx(128 * 118e-6)
    lsb = adc.lsb
    for k in (-10, -1, 0, 1, 10):
        code, sat = adc.convert(k * lsb + lsb / 2.0)
        assert int(code) == k
        assert not sat


def test_sar_flags_saturation_outside_its_range():
    adc = SarAdc(FdvpdParams(comparator_noise=0.0), _rng())
    for v, expect in ((0.0, False), (+1.0, True), (-1.0, True)):
        _, sat = adc.convert(v)
        assert bool(sat) is expect


def test_sar_output_is_clipped_not_wrapped():
    adc = SarAdc(FdvpdParams(comparator_noise=0.0), _rng())
    hi, _ = adc.convert(+1.0)
    lo, _ = adc.convert(-1.0)
    assert int(hi) == adc.n_codes // 2 - 1
    assert int(lo) == -adc.n_codes // 2


def test_comparator_noise_dithers_the_decision():
    adc = SarAdc(FdvpdParams(comparator_noise=100e-6), _rng())
    z = _rng(1).standard_normal(20000)
    codes, _ = adc.convert(np.zeros(20000) + adc.lsb / 2.0, z)
    # 100 uV rms against a 118 uV LSB -> the code must move around
    assert np.std(codes) == pytest.approx(100e-6 / adc.lsb, rel=0.15)
    assert np.mean(codes) == pytest.approx(0.0, abs=0.05)


def test_cdac_mismatch_keeps_the_quantiser_memoryless():
    adc = SarAdc(FdvpdParams(adc_sigma_cap=0.01, comparator_noise=0.0), _rng(2))
    assert np.all(np.diff(adc.levels) > 0)
    v = np.linspace(-5e-3, 5e-3, 101)
    a, _ = adc.convert(v)
    b, _ = adc.convert(v[::-1])
    assert np.array_equal(a, b[::-1])          # no hysteresis


def test_gross_cdac_mismatch_is_rejected():
    with pytest.raises(ValueError, match="non-monotonic"):
        SarAdc(FdvpdParams(adc_sigma_cap=2.0), _rng(0))


def test_adc_gain_trim_scales_the_lsb():
    adc = SarAdc(FdvpdParams(), _rng(), gain_trim=1.1)
    assert adc.lsb == pytest.approx(118e-6 * 1.1)


# ----------------------------------------------------- assembled FDVPD --
def _quiet_design(**kw) -> DesignParams:
    d = DesignParams()
    return replace(d, fdvpd=replace(d.fdvpd, comparator_noise=0.0, gamma=0.0,
                                    ktc_enabled=False, **kw))


def test_fdvpd_reads_zero_at_the_lock_point():
    """dt = T_frac*T_PD is the definition of lock: the residue must vanish."""
    d = _quiet_design()
    pd = Fdvpd(d, _rng())
    pd.prime(16)
    for t_frac in (0.05, 0.25, 0.5, 0.75, 0.95):
        s = pd.step(t_frac * d.t_pd, t_frac)
        assert abs(s.code) <= 1
        assert not s.saturated


def test_fdvpd_gain_is_one_lsb_per_lsb_of_time():
    """A timing error of one ADC LSB must move the code by exactly one."""
    d = _quiet_design()
    pd = Fdvpd(d, _rng())
    pd.prime(64)
    t_frac = 0.5
    lsb_t = d.fdvpd.adc_lsb_time
    # half-LSB offset so no sample lands exactly on a decision level, where
    # the tie would otherwise be broken by floating-point noise
    codes = [pd.step(t_frac * d.t_pd + (k + 0.5) * lsb_t, t_frac).code
             for k in range(-5, 6)]
    assert np.allclose(np.diff(codes), -1)     # a later edge = a negative code


def test_fdvpd_error_time_inverts_the_gain():
    d = _quiet_design()
    pd = Fdvpd(d, _rng())
    pd.prime(8)
    err = 3 * d.fdvpd.adc_lsb_time
    s = pd.step(0.5 * d.t_pd + err, 0.5)
    assert s.error_time == pytest.approx(-err, rel=0.02)


def test_fdvpd_saturates_outside_the_seven_bit_window():
    """The bang-bang regime the FLL exists to get out of (Sec. III-A step 5)."""
    d = _quiet_design()
    pd = Fdvpd(d, _rng())
    pd.prime(8)
    window = d.fdvpd.adc_full_scale / d.fdvpd.slew_rate    # ~18.9 ps
    assert not pd.step(0.5 * d.t_pd + 0.3 * window, 0.5).saturated
    assert pd.step(0.5 * d.t_pd + 2.0 * window, 0.5).saturated


def test_dac_code_is_clamped_to_the_array():
    d = _quiet_design()
    pd = Fdvpd(d, _rng())
    assert pd.dac_code(-0.5) == 0
    assert pd.dac_code(1.5) == pd.dac.n_codes - 1
    assert pd.dac_code(0.5) == pd.dac.n_codes // 2


def test_gain_mismatch_shifts_the_lock_point_with_the_fractional_word():
    """The dominant fractional-spur mechanism that PdGainCalibration removes."""
    d = _quiet_design(dac_gain_err=0.01)
    pd = Fdvpd(d, _rng())
    pd.prime(16)
    errs = []
    for t_frac in (0.1, 0.9):
        s = pd.step(t_frac * d.t_pd, t_frac)
        errs.append(s.error_time)
    # a gain error turns the fractional sawtooth into a proportional time error
    assert abs(errs[1]) > 5 * abs(errs[0])
    assert np.sign(errs[1]) == np.sign(errs[0])


def test_settling_makes_the_dac_code_dependent_on_its_predecessor():
    """Sec. III-E1: incomplete settling is a memory term, hence a spur."""
    d = _quiet_design(dac_settle_tau=2e-9)
    pd = Fdvpd(d, _rng())
    pd.prime(16)
    pd.step(0.9 * d.t_pd, 0.9)            # leave the DAC charged high
    after_high = pd.step(0.1 * d.t_pd, 0.1).code
    pd2 = Fdvpd(d, _rng())
    pd2.prime(16)
    pd2.step(0.1 * d.t_pd, 0.1)           # ...vs already settled low
    after_low = pd2.step(0.1 * d.t_pd, 0.1).code
    assert after_high != after_low


def test_prime_makes_the_noise_reproducible():
    d = default_design()
    out = []
    for _ in range(2):
        pd = Fdvpd(d, _rng(11))
        pd.prime(32)
        out.append([pd.step(0.5 * d.t_pd, 0.5).v_sampled for _ in range(32)])
    assert np.allclose(out[0], out[1])


def test_flicker_priming_needs_the_mean_window():
    """With a corner set, prime() must build the 1/f sequence."""
    d = replace(default_design(),
                fdvpd=replace(default_design().fdvpd, flicker_corner=520e3))
    pd = Fdvpd(d, _rng())
    pd.prime(4096, mean_t_on=100e-12)
    assert pd._v_flicker is not None
    assert len(pd._v_flicker) == 4096
    assert np.std(pd._v_flicker) > 0.0
