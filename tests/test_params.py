"""The derived design numbers must match the paper's text.

Every target below is quoted in Wu et al., JSSC 56(4):1254-1264, 2021; the
docstring of :mod:`fdvadpll.params` gives the section for each one.
"""

import math

import pytest

from fdvadpll import (DesignParams, default_design, fractional_design,
                      measured_fit_design)


# ---------------------------------------------------------------- top level --
def test_frequency_plan(design):
    """f_REF = 80 MHz, FCW = 42 -> 3.36 GHz, div-2 -> 1.68 GHz PD rate."""
    assert design.f_ref == pytest.approx(80e6)
    assert design.f_ckv == pytest.approx(3.36e9)
    assert design.f_pd == pytest.approx(1.68e9)
    assert design.t_pd == pytest.approx(595.24e-12, rel=1e-4)
    # the paper writes this channel as "FCW = 21 x 2"
    assert design.fcw_pd == pytest.approx(21.0)


def test_ramp_current_follows_from_sr_and_csar(design):
    """SR = 2*I_R/C_SAR with SR = 0.8 mV/ps and C_SAR = 550 fF -> 220 uA."""
    f = design.fdvpd
    assert f.slew_rate == pytest.approx(0.8e-3 / 1e-12)
    assert f.c_sar == pytest.approx(550e-15)
    assert f.i_ramp == pytest.approx(220e-6, rel=1e-9)
    assert f.slew_rate == pytest.approx(2.0 * f.i_ramp / f.c_sar)


# ------------------------------------------------------------------- I-DAC --
def test_dac_lsb(design):
    """10 b over one PD period: 581 fs = 465 uV."""
    assert design.dac_lsb_time == pytest.approx(581e-15, rel=2e-3)
    assert design.dac_lsb_volt == pytest.approx(465e-6, rel=2e-3)
    assert design.dac_lsb_volt == pytest.approx(
        design.dac_lsb_time * design.fdvpd.slew_rate)


def test_dac_full_scale_fits_the_ramp_range(design):
    """476 mV differential = 238 mV single ended, inside the 300 mV of Sec. III-B."""
    assert design.dac_full_scale_volt == pytest.approx(476e-3, rel=2e-3)
    single_ended = design.dac_full_scale_volt / 2.0
    assert single_ended == pytest.approx(238e-3, rel=2e-3)
    assert single_ended < 300e-3


def test_dac_segmentation(design):
    f = design.fdvpd
    assert f.dac_bits == 10
    assert f.dac_msb_thermo_bits == 4
    assert f.dac_bits - f.dac_msb_thermo_bits == 6


# --------------------------------------------------------------- SAR ADC --
def test_adc_lsb_and_range(design):
    """118 uV LSB = 147.5 fs; 7 b range = 15.1 mV, 'slightly more than 10 mV'."""
    f = design.fdvpd
    assert f.adc_lsb == pytest.approx(118e-6)
    assert f.adc_lsb_time == pytest.approx(147.5e-15, rel=1e-3)
    assert f.adc_full_scale == pytest.approx(15.1e-3, rel=2e-3)
    assert f.adc_full_scale > 10e-3


def test_subranging_gives_twelve_bits(design):
    """T_PD / 147.5 fs ~ 4035 ~ 12.0 bit -- the claim of Sec. II-B4."""
    assert design.effective_bits == pytest.approx(12.0, abs=0.05)
    assert design.t_pd / design.fdvpd.adc_lsb_time == pytest.approx(4035, rel=1e-3)


def test_subranging_beats_either_converter_alone(design):
    """The point of the architecture: 10 b DAC + 7 b ADC buy ~12 b, not 10 b."""
    f = design.fdvpd
    assert design.effective_bits > f.dac_bits
    assert design.effective_bits > f.adc_bits
    # ...but not the naive 10 + 7, because the ADC range overlaps the DAC LSB
    assert design.effective_bits < f.dac_bits + f.adc_bits


def test_adc_range_covers_more_than_one_dac_lsb(design):
    """Sub-ranging only works if the residue always lands inside the ADC."""
    assert design.fdvpd.adc_full_scale > design.dac_lsb_volt


# ------------------------------------------------------------- oscillator --
def test_period_jitter_round_trip(design):
    """sigma_dT is defined by L(1MHz) = f0^3 sigma^2 / (1MHz)^2."""
    sigma = design.dco.period_jitter_white()
    L = design.dco.f_center ** 3 * sigma ** 2 / (1e6 ** 2)
    assert 10.0 * math.log10(L) == pytest.approx(design.dco.pn_1mhz, abs=1e-9)


def test_reference_phase_noise_shape(design):
    """Flat floor well above the corner, +30 dB/decade below it."""
    far = design.ref.phase_noise(1e6)
    assert float(far) == pytest.approx(design.ref.pn_floor, abs=0.01)
    lo, hi = design.ref.phase_noise([10.0, 100.0])
    assert float(lo - hi) == pytest.approx(30.0, abs=0.2)


# ------------------------------------------------------------ loop filter --
def test_alpha_rho_map_to_bandwidth_and_damping(design):
    alpha, rho = design.loop.alpha_rho(design.f_ref)
    assert alpha * design.f_ref == pytest.approx(
        2.0 * math.pi * design.loop.bandwidth)
    assert alpha / (2.0 * math.sqrt(rho)) == pytest.approx(design.loop.damping)


def test_gear_shift_is_wider_than_the_tracking_bandwidth(design):
    assert design.loop.gear_bandwidth > design.loop.bandwidth


# ------------------------------------------------------------ design points --
def test_default_design_is_pure_physics():
    """No fitted term: the published equations only."""
    f = default_design().fdvpd
    assert f.flicker_corner == 0.0
    assert f.excess_noise_db == 0.0


def test_measured_fit_adds_exactly_one_parameter():
    """Sec. V blames the 1-100 kHz gap on ramp-generator flicker, nothing else."""
    d, m = default_design(), measured_fit_design()
    assert m.fdvpd.flicker_corner == pytest.approx(520e3)
    assert m.fdvpd.excess_noise_db == 0.0
    # everything else is untouched
    assert m.fdvpd.c_sar == d.fdvpd.c_sar
    assert m.fdvpd.comparator_noise == d.fdvpd.comparator_noise
    assert m.f_ckv == d.f_ckv


@pytest.mark.parametrize("bit,fcw_pd", [(2, 20.75), (5, 20.96875),
                                        (9, 20.998046875)])
def test_fractional_design_channels(bit, fcw_pd):
    """Only bit `frac_bit` asserted; Fig. 11 uses bit 5 -> 20.9688 x 2."""
    d = fractional_design(bit)
    assert d.fcw_pd == pytest.approx(fcw_pd, rel=1e-12)
    assert d.fcw == pytest.approx(2.0 * fcw_pd, rel=1e-12)


def test_with_fcw_leaves_the_original_alone(design):
    other = design.with_fcw(40.0)
    assert other.fcw == 40.0
    assert design.fcw == 42.0


def test_summary_renders(design):
    text = design.summary()
    assert "f_CKV" in text and "3.36000 GHz" in text
    assert "effective PD res." in text


def test_custom_design_stays_self_consistent():
    """Derived quantities must track a changed design point, not the defaults."""
    d = DesignParams(f_ref=100e6, fcw=40.0, fb_div=4)
    assert d.f_ckv == pytest.approx(4e9)
    assert d.f_pd == pytest.approx(1e9)
    assert d.fcw_pd == pytest.approx(10.0)
    assert d.dac_lsb_time == pytest.approx(d.t_pd / 1024)
