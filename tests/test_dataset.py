"""The design-space dataset generator.

Checks that the sampler covers what it claims, that feasibility is decided by
the constraints it documents, and that the targets are self-consistent -- not
that any particular design is good.
"""

import math
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "scripts"))

import make_dataset as md
from fdvadpll import noise


@pytest.fixture(scope="module")
def sample():
    return md.sample_space(n_pow2=7, seed=0)      # 128 points


# ---------------------------------------------------------------- sampling --
def test_sample_shape(sample):
    assert set(sample) == set(md.DIMS)
    for name, col in sample.items():
        assert len(col) == 128, name


def test_continuous_parameters_stay_inside_their_range(sample):
    for name, (lo, hi, _) in md.CONTINUOUS.items():
        col = sample[name]
        assert col.min() >= lo - 1e-9
        assert col.max() <= hi + 1e-9


def test_discrete_parameters_take_only_declared_values(sample):
    for name, values in md.DISCRETE.items():
        assert set(np.unique(sample[name])) <= set(values)


def test_log_scaled_parameters_are_log_uniform(sample):
    """Otherwise a decade-wide range would be dominated by its top decade."""
    col = np.sort(sample["c_sar_f"])              # declared 'log'
    lo, hi, how = md.CONTINUOUS["c_sar_f"]
    assert how == "log"
    median = float(np.median(col))
    assert median == pytest.approx(math.sqrt(lo * hi), rel=0.25)


def test_frac_bit_encoding(sample):
    frac = sample["frac_bit"]
    assert set(np.unique(frac)) <= {-1} | set(md.FRAC_BITS)
    share_integer = float((frac < 0).mean())
    assert share_integer == pytest.approx(md.P_INTEGER_N, abs=0.08)


def test_sampling_is_reproducible():
    a = md.sample_space(6, seed=3)
    b = md.sample_space(6, seed=3)
    for k in a:
        assert np.array_equal(a[k], b[k])


def test_seed_changes_the_sample():
    a = md.sample_space(6, seed=3)
    b = md.sample_space(6, seed=4)
    assert not np.array_equal(a["c_sar_f"], b["c_sar_f"])


def test_sobol_covers_each_dimension_evenly(sample):
    """The reason for a low-discrepancy sequence rather than uniform random."""
    lo, hi, _ = md.CONTINUOUS["loop_damping"]
    u = (sample["loop_damping"] - lo) / (hi - lo)
    counts, _ = np.histogram(u, bins=8, range=(0, 1))
    assert counts.min() >= 128 // 8 - 2          # near-perfect stratification


# ------------------------------------------------------------ the design --
def test_build_design_round_trips_the_parameters(sample):
    row = {k: v[0] for k, v in sample.items()}
    d = md.build_design(row)
    assert d.f_ref == md.F_REF
    assert d.fcw == pytest.approx(row["fcw"])
    assert d.fb_div == int(row["fb_div"])
    assert d.fdvpd.c_sar == pytest.approx(row["c_sar_f"])
    assert d.fdvpd.adc_bits == int(row["adc_bits"])
    assert d.loop.bandwidth == pytest.approx(row["loop_bw_hz"])


def test_oscillator_is_centred_on_the_channel(sample):
    """Regression: f_center used to stay at the prototype's 3.36 GHz.

    Sweeping fcw over 2-8 GHz while the oscillator sat at 3.36 GHz handed every
    design a multi-gigahertz acquisition problem unrelated to the parameter
    under study.
    """
    for i in range(16):
        row = {k: v[i] for k, v in sample.items()}
        d = md.build_design(row)
        assert d.dco.f_center == pytest.approx(d.f_ckv)


# ------------------------------------------------------------ constraints --
def test_check_passes_the_prototype():
    from fdvadpll import default_design
    d = default_design()
    margin = md.capture_margin(d)["pd_capture_margin"]
    assert md.check(d, margin, "wide") == ""
    assert md.check(d, margin, "realisable") == ""


def test_check_catches_a_broken_subranging_design():
    from dataclasses import replace
    from fdvadpll import default_design
    d = default_design()
    # a 12 b DAC with a 5 b ADC: the residue no longer fits the converter
    bad = replace(d, fdvpd=replace(d.fdvpd, dac_bits=8, adc_bits=5,
                                   adc_lsb=30e-6))
    assert "sub-ranging" in md.check(bad, 100.0, "wide")


def test_check_catches_an_oversized_ramp_swing():
    from dataclasses import replace
    from fdvadpll import default_design
    d = default_design()
    # c_sar trimmed so the ramp current stays realisable and this is the only
    # constraint left to trip
    bad = replace(d, fdvpd=replace(d.fdvpd, slew_rate=20e9, c_sar=200e-15))
    assert bad.fdvpd.i_ramp < md.I_RAMP_MAX
    assert "ramp swing" in md.check(bad, 100.0, "wide")


def test_realisable_profile_is_stricter_than_wide():
    from dataclasses import replace
    from fdvadpll import default_design
    d = default_design()
    wide_bw = replace(d, loop=replace(d.loop,
                                      bandwidth=md.MAX_REALISABLE_BW * 2))
    assert md.check(wide_bw, 100.0, "wide") == ""
    assert "locking envelope" in md.check(wide_bw, 100.0, "realisable")
    assert md.check(d, 1.0, "wide") == ""
    assert "rail" in md.check(d, 1.0, "realisable")


# -------------------------------------------------------- capture margin --
def test_capture_margin_on_the_prototype():
    from fdvadpll import measured_fit_design
    m = md.capture_margin(measured_fit_design())
    d = measured_fit_design()
    assert m["pd_capture_s"] == pytest.approx(
        d.fdvpd.adc_full_scale / d.fdvpd.slew_rate)
    assert m["pd_residual_s"] > 0
    # the prototype is comfortably inside its own detector
    assert m["pd_capture_margin"] > md.MIN_CAPTURE_MARGIN


def test_a_wider_loop_leaves_less_residual_for_the_detector():
    """More loop gain corrects more of the oscillator, so the margin grows."""
    from dataclasses import replace
    from fdvadpll import measured_fit_design
    d = measured_fit_design()
    narrow = replace(d, loop=replace(d.loop, bandwidth=100e3))
    wide = replace(d, loop=replace(d.loop, bandwidth=1e6))
    assert (md.capture_margin(wide)["pd_residual_s"]
            < md.capture_margin(narrow)["pd_residual_s"])


def test_a_noisier_oscillator_eats_the_margin():
    from dataclasses import replace
    from fdvadpll import measured_fit_design
    d = measured_fit_design()
    loud = replace(d, dco=replace(d.dco, pn_1mhz=-115.0))
    assert (md.capture_margin(loud)["pd_capture_margin"]
            < md.capture_margin(d)["pd_capture_margin"])


# ----------------------------------------------------------- the targets --
def test_evaluate_is_self_consistent():
    from fdvadpll import measured_fit_design
    d = measured_fit_design()
    out = md.evaluate(d, frac_bit=5)

    budget = noise.pd_noise_floor(d, 5)
    assert out["L_pd_total_dbc"] == pytest.approx(
        10 * math.log10(budget.pd_total))
    assert out["fom_db"] == pytest.approx(
        noise.fom(out["jitter_s"], out["p_total_w"]))
    assert out["p_total_w"] == pytest.approx(out["p_fdvpd_w"] + md.P_REST_W)
    assert out["dominant_pd_term"] in ("ramp", "ktc", "comparator",
                                       "quantisation")


def test_evaluate_matches_the_published_numbers():
    """The generator must reproduce the paper through the same code path."""
    from fdvadpll import measured_fit_design
    out = md.evaluate(measured_fit_design(), frac_bit=None)
    assert out["jitter_s"] * 1e15 == pytest.approx(82.0, abs=4.0)
    assert out["L_out_100khz_dbc"] == pytest.approx(-116.0, abs=2.0)


def test_probe_columns_lie_on_the_profile():
    from fdvadpll import measured_fit_design
    d = measured_fit_design()
    out = md.evaluate(d, None)
    total = noise.output_phase_noise(d, md.F_GRID, None)
    for probe, key in ((10e3, "L_out_10khz_dbc"), (1e6, "L_out_1mhz_dbc")):
        want = 10 * math.log10(float(np.interp(probe, md.F_GRID, total)))
        assert out[key] == pytest.approx(want)


def test_derived_columns_agree_with_the_design():
    from fdvadpll import measured_fit_design
    d = measured_fit_design()
    got = md.derived(d)
    assert got["f_ckv_hz"] == pytest.approx(d.f_ckv)
    assert got["i_ramp_a"] == pytest.approx(d.fdvpd.i_ramp)
    assert got["effective_bits"] == pytest.approx(d.effective_bits)
    assert got["adc_fs_v"] == pytest.approx(d.fdvpd.adc_full_scale)


@pytest.mark.slow
def test_cross_check_reports_a_locked_prototype():
    from fdvadpll import measured_fit_design
    out = md.cross_check(measured_fit_design(), seed=1, n_cycles=1 << 15)
    assert out["sim_locked"] is True
    assert out["sim_cycle_slips"] == 0
    assert out["sim_jitter_s"] * 1e15 == pytest.approx(82.0, abs=15.0)
