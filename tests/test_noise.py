"""The s-domain model of Sec. II-B, equations (1)-(13).

The headline test is :func:`test_measured_fit_reproduces_the_published_numbers`:
with the single fitted flicker corner the model must land on the jitter and
in-band noise the paper actually measured.
"""

import math

import numpy as np
import pytest

from fdvadpll import noise
from fdvadpll.params import KT, DesignParams, default_design, measured_fit_design


F_GRID = np.logspace(3, math.log10(40e6), 4000)


# ------------------------------------------- equation (1): gating window --
def test_integer_channel_duty_cycle_is_the_fixed_on_time(design):
    """Integer-N: the window collapses to the 100 ps V_OS margin."""
    duty = noise.gating_duty_cycle(design, None, t_on_integer=100e-12)
    assert duty == pytest.approx(100e-12 * design.f_ref)


@pytest.mark.parametrize("bit", [2, 5, 9, 16])
def test_fractional_duty_cycle_formula(design, bit):
    """g_bar = (1/2 - 2^-(Nq+1)) * T_PD, averaged over the 2^Nq-cycle pattern."""
    duty = noise.gating_duty_cycle(design, bit)
    expected = (0.5 - 0.5 ** (bit + 1)) * design.t_pd * design.f_ref
    assert duty == pytest.approx(expected)
    # explicit average over the sawtooth the accumulator actually walks through
    widths = np.arange(2 ** bit) / 2 ** bit * design.t_pd
    assert duty == pytest.approx(widths.mean() * design.f_ref)


def test_duty_cycle_doubles_without_the_div2(design):
    """The div-2 doubles the PD period, hence the window -- the +3 dB of Sec. II-B."""
    with_div = noise.gating_duty_cycle(design, 5, include_div2=True)
    without = noise.gating_duty_cycle(design, 5, include_div2=False)
    assert with_div / without == pytest.approx(design.fb_div)
    assert 10 * math.log10(with_div / without) == pytest.approx(3.01, abs=0.01)


@pytest.mark.parametrize("bit", [3, 5])
def test_gating_harmonics_satisfies_parseval(bit):
    """sum_n |G_n|^2 must converge to the mean of g(t)^2  --  eqns (3)-(5)."""
    exact = noise.gating_harmonics(bit)
    assert exact == pytest.approx((0.5 - 0.5 ** (bit + 1)))
    truncated = noise.gating_harmonics(bit, n_harmonics=3000)
    assert truncated == pytest.approx(exact, rel=0.01)
    # and it must approach from below as more harmonics are kept
    coarse = noise.gating_harmonics(bit, n_harmonics=50)
    assert coarse < truncated <= exact * 1.001


# --------------------------------------- equations (2), (6), (10): gains --
def test_ramp_current_psd_is_4kTgamma_2gm(design):
    f = design.fdvpd
    assert noise.ramp_current_psd(design) == pytest.approx(
        4.0 * KT * f.gamma * 2.0 * f.gm_ramp)


def test_k_ramp_and_k_slope(design):
    f = design.fdvpd
    assert noise.k_ramp(design) == pytest.approx(
        f.i_ramp / (math.pi * design.n_mult))
    assert noise.k_slope(design) == pytest.approx(
        2.0 * f.i_ramp / (2.0 * math.pi * design.f_ckv * f.c_sar))
    # K_SR is just the slew rate referred to output phase
    assert noise.k_slope(design) == pytest.approx(
        f.slew_rate / (2.0 * math.pi * design.f_ckv))


def test_gains_are_independent_of_the_feedback_divider(design):
    """A timing error is a timing error: K_R and K_SR do not see the div-2."""
    other = DesignParams(fb_div=4, fcw=design.fcw)
    assert noise.k_slope(other) == pytest.approx(noise.k_slope(design))
    assert noise.k_ramp(other) == pytest.approx(noise.k_ramp(design))


# ------------------------------- equations (9)-(13): the four noise terms --
def test_ktc_closed_form(design):
    f = design.fdvpd
    assert noise.L_ktc(design) == pytest.approx(
        2.0 * math.pi ** 2 * design.f_ckv ** 2 / design.f_ref
        * KT * f.c_sar / f.i_ramp ** 2)


def test_quantisation_closed_form(design):
    """pi^2/3 * f_CKV^2/f_REF * LSB_t^2, written with the time LSB."""
    lsb = design.fdvpd.adc_lsb_time
    assert noise.L_quantisation(design) == pytest.approx(
        math.pi ** 2 / 3.0 * design.f_ckv ** 2 / design.f_ref * lsb ** 2)


def test_quantisation_drops_6db_per_bit(design):
    """One extra bit halves the time LSB."""
    lsb = design.fdvpd.adc_lsb_time
    a = noise.L_quantisation(design, lsb)
    b = noise.L_quantisation(design, lsb / 2.0)
    assert 10 * math.log10(a / b) == pytest.approx(6.02, abs=0.01)


def test_twelve_bit_over_t_pd_equals_eleven_bit_over_t_ckv(design):
    """The div-2 bookkeeping note in the module docstring."""
    a = noise.L_quantisation(design, design.t_pd / 2 ** 12)
    b = noise.L_quantisation(design, design.t_ckv / 2 ** 11)
    assert a == pytest.approx(b)


def test_comparator_noise_scales_with_v_tn_squared(design):
    from dataclasses import replace
    base = noise.L_comparator(design)
    quiet = noise.L_comparator(
        replace(design, fdvpd=replace(design.fdvpd, comparator_noise=50e-6)))
    assert base / quiet == pytest.approx(4.0)


# ------------------------------------------------------- the full budget --
def test_budget_totals_its_parts(design):
    b = noise.pd_noise_floor(design, 5)
    assert b.pd_total == pytest.approx(b.ramp + b.ktc + b.comparator
                                       + b.quantisation)
    d = b.as_dbc()
    assert d["pd_total"] == pytest.approx(10 * math.log10(b.pd_total))


def test_quantisation_is_not_the_limit(design):
    """Sub-ranging was worth it: the ADC LSB sits well below the thermal floor."""
    b = noise.pd_noise_floor(design, 5)
    assert b.quantisation < 0.1 * b.pd_total


def test_fractional_channel_is_noisier_than_integer(design):
    """A wider ramp window integrates more I_R noise  --  equation (1)."""
    integer = noise.pd_noise_floor(design, None).pd_total
    frac = noise.pd_noise_floor(design, 5).pd_total
    assert frac > integer


def test_only_the_ramp_term_cares_about_the_fractional_word(design):
    a = noise.pd_noise_floor(design, None)
    b = noise.pd_noise_floor(design, 5)
    assert a.ktc == pytest.approx(b.ktc)
    assert a.comparator == pytest.approx(b.comparator)
    assert a.quantisation == pytest.approx(b.quantisation)
    assert b.ramp > a.ramp


def test_excess_noise_scales_thermal_but_not_quantisation(design):
    from dataclasses import replace
    loud = replace(design, fdvpd=replace(design.fdvpd, excess_noise_db=3.0))
    a, b = noise.pd_noise_floor(design, 5), noise.pd_noise_floor(loud, 5)
    assert b.ramp / a.ramp == pytest.approx(10 ** 0.3)
    assert b.ktc / a.ktc == pytest.approx(10 ** 0.3)
    assert b.quantisation == pytest.approx(a.quantisation)


def test_fdvpd_beats_an_equivalent_charge_pump_by_more_than_10db(design):
    """The >10 dB claim of Sec. II-B1, equation (8) vs equation (7)."""
    cp = noise.L_charge_pump_reference(design)
    fdv = noise.L_ramp_current(design, 5)
    assert 10 * math.log10(cp / fdv) > 10.0


def test_flicker_shaping_crosses_the_thermal_floor_at_the_corner():
    """PD noise rises as (1 + f_c/f): exactly 3 dB up at f = f_c."""
    d = measured_fit_design()
    fc = d.fdvpd.flicker_corner
    _, total = noise.pd_noise_profile(d, np.array([fc, 100 * fc]), 5)
    flat = noise.pd_noise_floor(d, 5).pd_total
    assert total[0] / flat == pytest.approx(2.0, rel=0.02)
    assert total[1] / flat == pytest.approx(1.0, rel=0.02)


def test_profile_parts_sum_to_the_total(design):
    parts, total = noise.pd_noise_profile(design, F_GRID, 5)
    assert np.allclose(sum(parts.values()), total)


# ------------------------------------------------- closed-loop behaviour --
def test_loop_transfer_is_complementary(design):
    h_lp, h_hp = noise.loop_transfer(design, F_GRID)
    assert np.allclose(h_lp + h_hp, 1.0)


def test_loop_transfer_limits(design):
    h_lp, h_hp = noise.loop_transfer(design, np.array([1.0, 1e9]))
    assert abs(h_lp[0]) == pytest.approx(1.0, rel=1e-6)   # DC: tracks
    assert abs(h_hp[0]) < 1e-9                            # DC: rejects
    assert abs(h_hp[1]) == pytest.approx(1.0, rel=1e-3)   # HF: passes DCO
    assert abs(h_lp[1]) < 1e-2


def test_unity_gain_lands_at_the_requested_bandwidth(design):
    """|H_lp| = |H_hp| where the loop gain is unity."""
    f = np.logspace(4, 7, 20001)
    h_lp, h_hp = noise.loop_transfer(design, f)
    cross = f[np.argmin(np.abs(np.abs(h_lp) - np.abs(h_hp)))]
    assert cross == pytest.approx(design.loop.bandwidth, rel=0.05)


def test_output_breakdown_sums_to_the_total(design):
    total, parts = noise.output_phase_noise(design, F_GRID, 5, breakdown=True)
    shaped = {k: v for k, v in parts.items() if k != "DCO (free running)"}
    assert np.allclose(sum(shaped.values()), total)


def test_loop_suppresses_the_dco_in_band(design):
    """Deep in band the oscillator is corrected away; the PD is what is left."""
    f = np.array([10e3])
    total, parts = noise.output_phase_noise(design, f, 5, breakdown=True)
    free = parts["DCO (free running)"][0]
    assert parts["DCO"][0] < 1e-3 * free          # >30 dB of rejection
    assert parts["DCO"][0] < 0.1 * total[0]
    pd = sum(parts[k][0] for k in ("ramp current", "kT/C", "comparator",
                                   "quantisation"))
    assert pd > 0.5 * total[0]


def test_dco_takes_over_around_the_loop_bandwidth(design):
    """The DCO share grows from negligible in band to the largest single term."""
    f = np.array([10e3, 1e6])
    _, parts = noise.output_phase_noise(design, f, 5, breakdown=True)
    share = parts["DCO"] / sum(v for k, v in parts.items()
                               if k != "DCO (free running)")
    assert share[1] > 10 * share[0]
    assert share[1] == max(
        float(v[1]) for k, v in parts.items() if k != "DCO (free running)") \
        / float(sum(v[1] for k, v in parts.items()
                    if k != "DCO (free running)"))


def test_type2_tail_leaves_pd_noise_unfiltered(design):
    """Why :class:`DigitalLoopFilter` offers extra IIR poles.

    A bare type-II loop rolls ``H_lp`` off at only -20 dB/decade, the same slope
    as the free-running DCO, so far outside the bandwidth the PD contribution
    never stops tracking the oscillator -- it just sits a fixed number of dB
    below it.  Real ADPLLs add poles to break that.
    """
    f = np.array([5e6, 4e7])
    _, parts = noise.output_phase_noise(design, f, 5, breakdown=True)
    ratio = parts["ramp current"] / parts["DCO"]
    assert ratio[1] == pytest.approx(ratio[0], rel=0.02)
    # and the composite therefore falls at 20 dB/decade out there
    total = noise.output_phase_noise(design, f, 5)
    decades = math.log10(f[1] / f[0])
    assert 10 * math.log10(total[0] / total[1]) / decades == pytest.approx(
        20.0, abs=0.5)


def test_reference_noise_is_multiplied_by_n_squared(design):
    f = np.array([1e4])
    _, parts = noise.output_phase_noise(design, f, 5, breakdown=True)
    raw = 10 ** (design.ref.phase_noise(f) / 10.0)
    h_lp, _ = noise.loop_transfer(design, f)
    assert parts["reference"][0] == pytest.approx(
        raw[0] * design.n_mult ** 2 * abs(h_lp[0]) ** 2)


def test_widening_the_loop_trades_pd_noise_for_dco_noise(design):
    """The classic optimum: there is a bandwidth that minimises jitter."""
    from dataclasses import replace
    jitters = []
    for bw in (100e3, 500e3, 3e6):
        d = replace(design, loop=replace(design.loop, bandwidth=bw))
        L = noise.output_phase_noise(d, F_GRID, 5)
        jitters.append(noise.integrated_jitter(F_GRID, L, d.f_ckv))
    assert jitters[1] < jitters[0] and jitters[1] < jitters[2]


# ------------------------------------------------------- jitter and FoM --
def test_integrated_jitter_matches_the_closed_form():
    f = np.logspace(3, 8, 40001)
    A, f0 = 1e-2, 3.36e9
    expected = (math.sqrt(2 * A * (1 / 10e3 - 1 / 40e6))
                / (2 * math.pi * f0))
    assert noise.integrated_jitter(f, A / f ** 2, f0) == pytest.approx(
        expected, rel=1e-3)


def test_integrated_jitter_rejects_a_degenerate_band():
    with pytest.raises(ValueError):
        noise.integrated_jitter(np.array([1e3, 1e4]), np.ones(2), 3.36e9,
                                f_lo=1e6, f_hi=2e6)


def test_fom_definition():
    assert noise.fom(100e-15, 9.2e-3) == pytest.approx(
        20 * math.log10(100e-15) + 10 * math.log10(9.2), abs=1e-9)


def test_fom_improves_6db_when_jitter_halves():
    assert noise.fom(50e-15, 9.2e-3) - noise.fom(100e-15, 9.2e-3) == \
        pytest.approx(-6.02, abs=0.01)


def test_normalise_spur_is_identity_at_the_reference_carrier():
    assert noise.normalise_spur(-60.0, 3.24e9) == pytest.approx(-60.0)
    assert noise.normalise_spur(-60.0, 6.48e9) == pytest.approx(-66.02, abs=0.01)


# ------------------------------------------------ the published data points --
def test_default_design_underpredicts_the_measurement(design):
    """Sec. V: without the ramp flicker term the model is optimistic in band."""
    L = noise.output_phase_noise(design, F_GRID, None)
    at_110k = 10 * math.log10(np.interp(110e3, F_GRID, L))
    assert at_110k < -116.0          # better than measured, i.e. optimistic
    assert at_110k > -125.0          # but not absurdly so


def test_measured_fit_reproduces_the_published_numbers():
    """The whole point of the fit: 82 fs / -116 dBc/Hz integer, 101 fs fractional.

    Targets are the measured values of Fig. 11 and Sec. V.
    """
    d = measured_fit_design()

    L_int = noise.output_phase_noise(d, F_GRID, None)
    j_int = noise.integrated_jitter(F_GRID, L_int, d.f_ckv)
    assert j_int * 1e15 == pytest.approx(82.0, abs=3.0)

    at_110k = 10 * math.log10(np.interp(110e3, F_GRID, L_int))
    assert at_110k == pytest.approx(-116.0, abs=1.0)

    L_frac = noise.output_phase_noise(d, F_GRID, 5)
    j_frac = noise.integrated_jitter(F_GRID, L_frac, d.f_ckv)
    assert j_frac * 1e15 == pytest.approx(101.0, abs=6.0)

    # ...and the FoM the paper reports for 9.2 mW
    assert noise.fom(j_int, 9.2e-3) == pytest.approx(-252.0, abs=1.5)


# ------------------------------------------------------- Fig. 5 power model --
def test_fdvpd_power_terms_add_up(design):
    f = design.fdvpd
    p = noise.fdvpd_power(design, f.i_ramp, f.comparator_noise, f.c_sar)
    assert p["total"] == pytest.approx(p["idac"] + p["cdyn"] + p["comparator"]
                                       + p["logic"])
    assert all(v > 0 for v in p.values())


def test_comparator_power_follows_the_energy_noise_product(design):
    """30 nJ*uV^2 per conversion: halving v_TN costs 4x the power."""
    f = design.fdvpd
    a = noise.fdvpd_power(design, f.i_ramp, 100e-6, f.c_sar)["comparator"]
    b = noise.fdvpd_power(design, f.i_ramp, 50e-6, f.c_sar)["comparator"]
    assert b / a == pytest.approx(4.0)


def test_noise_floor_falls_with_power(design):
    """Fig. 5: more FDVPD power buys a lower in-band floor, monotonically."""
    powers = np.array([0.5e-3, 1e-3, 2e-3, 4e-3, 8e-3])
    L, info = noise.noise_vs_power(design, 12, powers)
    assert np.all(np.isfinite(L))
    assert np.all(np.diff(L) < 0)
    assert np.all(np.diff(info["i_ramp"]) > 0)


def test_more_pd_bits_lower_the_floor(design):
    """...until the thermal terms take over, so the gain is less than 6 dB/bit."""
    powers = np.array([2e-3])
    l10 = noise.noise_vs_power(design, 10, powers)[0][0]
    l12 = noise.noise_vs_power(design, 12, powers)[0][0]
    assert l12 < l10
    assert l10 - l12 < 12.0


def test_infeasible_power_budget_is_flagged(design):
    L, info = noise.noise_vs_power(design, 12, np.array([1e-6]))
    assert not info["feasible"][0]
    assert math.isnan(L[0])
