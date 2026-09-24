"""Generate a design-space dataset from the analytic model.

Sobol-samples the FDVPD / loop / oscillator design space, evaluates the
s-domain model of Sec. II-B at every point, and writes one CSV row per design
with the parameters as features and the noise, jitter, power and FoM as
targets.  A provenance file records the ranges, the fixed choices and the
limitations, so the CSV is interpretable on its own.

    python scripts/make_dataset.py                     # 16384 points, ~1 min
    python scripts/make_dataset.py --n 65536 --sim 32  # bigger, + sim check

Why Sobol rather than a grid: a grid over 18 parameters is hopeless, and
uniform random sampling clumps.  A low-discrepancy sequence covers the space
evenly at any sample count, which is what a surrogate model wants.

Three things to know before using the output
--------------------------------------------
*Infeasible rows are kept, not dropped.*  A design whose ADC range cannot cover
one DAC LSB is a real and interesting part of the space -- it is where
sub-ranging stops working.  ``feasible`` flags those rows and ``reject_reason``
says which constraint failed, so you can filter, or learn the constraint.

*The labels come from the analytic model, not the simulator.*  That is what
makes 10**4 points affordable (milliseconds each, against ~6 s for an
event-driven run).  ``--sim`` cross-checks a random subset with the full
simulator and writes the result into the same rows, so the dataset carries its
own estimate of label accuracy.  On the default design point the two agree to
about 3 %.

*The oscillator's own power is not modelled.*  ``dco_pn_1mhz_dbc`` is a free
parameter here, but in silicon a quieter oscillator costs more power, and no
such relation exists in this package.  Total power is therefore the *computed*
FDVPD power plus a fixed allowance for everything else, taken from the
prototype's measured breakdown.  Treat ``fom_db`` as a detector figure of
merit, not a whole-synthesiser one.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import sys
import time
from dataclasses import replace

import numpy as np
from scipy.stats import qmc

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from fdvadpll import DesignParams, FdvPll, __version__, noise

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# --------------------------------------------------------------------------
#  Fixed choices
# --------------------------------------------------------------------------
#: The reference stays at the prototype's 80 MHz.  Sweeping it would make the
#: fixed 10 kHz - 40 MHz integration band mean different things in different
#: rows, and the jitter column would stop being comparable.
F_REF = 80e6

#: Integration band for the jitter column, same for every row.
F_LO, F_HI = 10e3, 40e6

#: Offsets the output profile is sampled at, as extra target columns.
PROBES = (10e3, 100e3, 1e6, 10e6)

#: Differential swing the ramp can produce.  The prototype uses 476 mV out of a
#: 300 mV single-ended (600 mV differential) range; allow a little more.
RAMP_HEADROOM_V = 1.2

#: Ramp currents outside this are not realisable at this slew rate.
I_RAMP_MIN, I_RAMP_MAX = 10e-6, 5e-3

#: Extra constraints applied by the 'realisable' profile.  Both are empirical,
#: measured by cross-checking Sobol samples against the event-driven simulator
#: -- neither follows from the s-domain model, which has no notion of lock.
MIN_CAPTURE_MARGIN = 6.0
MAX_REALISABLE_BW = 1.0e6

#: Everything that is not the FDVPD, from the prototype's Fig. 17 breakdown:
#: oscillator + divider + loop filter + the rest, held fixed.  See the module
#: docstring on why this is an allowance and not a model.
_P = DesignParams().power
P_REST_W = (_P.share_dco_div2 + _P.share_dlf + _P.share_rest) * _P.p_total_measured

#: Frequency grid the profile is integrated on.
F_GRID = np.logspace(3, math.log10(F_HI), 2000)


def _show(path: pathlib.Path) -> str:
    """Repo-relative if it is inside the repo, absolute otherwise."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


# --------------------------------------------------------------------------
#  The sampled space
# --------------------------------------------------------------------------
#: name -> (low, high, scale).  'log' samples log-uniformly, 'lin' uniformly.
CONTINUOUS = {
    "fcw":                 (25.0, 100.0, "lin"),    # f_CKV 2 .. 8 GHz
    "sr_v_per_s":          (0.2e9, 3.2e9, "log"),   # 0.2 .. 3.2 mV/ps
    "c_sar_f":             (100e-15, 3000e-15, "log"),
    "adc_lsb_v":           (30e-6, 400e-6, "log"),
    "comparator_noise_v":  (25e-6, 400e-6, "log"),
    "gamma":               (0.7, 2.0, "lin"),
    "gm_id":               (5.0, 20.0, "lin"),
    "flicker_corner_hz":   (10e3, 2e6, "log"),
    "dco_pn_1mhz_dbc":     (-140.0, -120.0, "lin"),
    "dco_fcorner_hz":      (10e3, 500e3, "log"),
    "ref_pn_floor_dbc":    (-170.0, -150.0, "lin"),
    "ref_fcorner_hz":      (100.0, 10e3, "log"),
    "loop_bw_hz":          (50e3, 4e6, "log"),
    "loop_damping":        (0.5, 2.0, "lin"),
}

#: name -> the values it may take
DISCRETE = {
    "fb_div":   (1, 2, 4),
    "dac_bits": (8, 9, 10, 11, 12),
    "adc_bits": (5, 6, 7, 8, 9),
}

#: probability that a row is an integer-N channel; otherwise a fractional bit
#: is drawn uniformly from FRAC_BITS
P_INTEGER_N = 0.25
FRAC_BITS = tuple(range(2, 17))

DIMS = list(CONTINUOUS) + list(DISCRETE) + ["frac_bit"]


def _scale(u: np.ndarray, lo: float, hi: float, how: str) -> np.ndarray:
    if how == "log":
        return np.exp(np.log(lo) + u * (math.log(hi) - math.log(lo)))
    return lo + u * (hi - lo)


def sample_space(n_pow2: int, seed: int) -> dict[str, np.ndarray]:
    """Sobol sample of the design space.  Returns a column-per-parameter dict."""
    sob = qmc.Sobol(d=len(DIMS), scramble=True, seed=seed)
    u = sob.random_base2(n_pow2)

    cols: dict[str, np.ndarray] = {}
    for i, (name, (lo, hi, how)) in enumerate(CONTINUOUS.items()):
        cols[name] = _scale(u[:, i], lo, hi, how)

    j = len(CONTINUOUS)
    for name, values in DISCRETE.items():
        idx = np.minimum((u[:, j] * len(values)).astype(int), len(values) - 1)
        cols[name] = np.asarray(values, dtype=int)[idx]
        j += 1

    # integer-N is encoded as frac_bit = -1
    uf = u[:, j]
    frac = np.full(len(uf), -1, dtype=int)
    is_frac = uf >= P_INTEGER_N
    pick = ((uf[is_frac] - P_INTEGER_N) / (1.0 - P_INTEGER_N) * len(FRAC_BITS))
    frac[is_frac] = np.asarray(FRAC_BITS, dtype=int)[
        np.minimum(pick.astype(int), len(FRAC_BITS) - 1)]
    cols["frac_bit"] = frac
    return cols


# --------------------------------------------------------------------------
#  One design point
# --------------------------------------------------------------------------
def build_design(row: dict) -> DesignParams:
    d = DesignParams()
    return replace(
        d,
        f_ref=F_REF,
        fcw=float(row["fcw"]),
        fb_div=int(row["fb_div"]),
        fdvpd=replace(
            d.fdvpd,
            slew_rate=float(row["sr_v_per_s"]),
            c_sar=float(row["c_sar_f"]),
            dac_bits=int(row["dac_bits"]),
            adc_bits=int(row["adc_bits"]),
            adc_lsb=float(row["adc_lsb_v"]),
            comparator_noise=float(row["comparator_noise_v"]),
            gamma=float(row["gamma"]),
            gm_id=float(row["gm_id"]),
            flicker_corner=float(row["flicker_corner_hz"]),
        ),
        # The oscillator is centred on the channel being synthesised.  Leaving
        # f_center at the prototype's 3.36 GHz while fcw moves over 2-8 GHz
        # would hand every design a multi-gigahertz acquisition problem that
        # has nothing to do with the parameter being studied.
        dco=replace(d.dco, f_center=float(row["fcw"]) * F_REF,
                    pn_1mhz=float(row["dco_pn_1mhz_dbc"]),
                    f_corner_flicker=float(row["dco_fcorner_hz"])),
        ref=replace(d.ref, f_ref=F_REF, pn_floor=float(row["ref_pn_floor_dbc"]),
                    f_corner=float(row["ref_fcorner_hz"])),
        loop=replace(d.loop, bandwidth=float(row["loop_bw_hz"]),
                     damping=float(row["loop_damping"])),
    )


def capture_margin(design: DesignParams) -> dict:
    """How much room the detector has before it rails.

    The phase detector only captures ``adc_full_scale/SR`` of timing error.
    What it has to hold is the residual the loop does *not* correct: the
    oscillator and reference noise high-pass filtered by the loop.  The ratio
    of the two is the headroom, quoted here as capture window over six sigma.

    This is the constraint the analytic model is blind to.  ``noise.py`` will
    happily report a jitter figure for a design whose detector would be railed
    from the first cycle, because nothing in the s-domain analysis knows the
    ADC has ends.  Below a margin of about 5 the simulator stops locking.
    """
    h_lp, h_hp = noise.loop_transfer(design, F_GRID)
    residual = ((noise.dco_free_running_pn(design, F_GRID)
                 + 10.0 ** (design.ref.phase_noise(F_GRID) / 10.0)
                 * design.n_mult ** 2)
                * np.abs(h_hp) ** 2)
    sigma_t = noise.integrated_jitter(F_GRID, residual, design.f_ckv,
                                      1e3, F_HI)
    capture = design.fdvpd.adc_full_scale / design.fdvpd.slew_rate
    return {
        "pd_capture_s": capture,
        "pd_residual_s": sigma_t,
        "pd_capture_margin": capture / (6.0 * sigma_t) if sigma_t > 0 else
                             float("inf"),
    }


def check(design: DesignParams, margin: float, profile: str) -> str:
    """Return the first constraint the design violates, or '' if it is fine."""
    f = design.fdvpd
    if f.i_ramp < I_RAMP_MIN or f.i_ramp > I_RAMP_MAX:
        return "ramp current outside a realisable range"
    if design.dac_full_scale_volt > RAMP_HEADROOM_V:
        return "DAC full scale exceeds the ramp swing"
    if f.adc_full_scale < 2.0 * design.dac_lsb_volt:
        return "ADC range cannot cover one DAC LSB (sub-ranging broken)"
    if design.loop.bandwidth > design.f_ref / 10.0:
        return "loop bandwidth above f_REF/10"
    if profile == "realisable":
        # Empirical, from cross-checking Sobol samples against the simulator:
        # below this margin the detector rails, and above ~1 MHz of loop
        # bandwidth designs stop locking for reasons the s-domain model does
        # not capture.  See the module docstring.
        if margin < MIN_CAPTURE_MARGIN:
            return "phase detector would rail on the uncorrected residual"
        if design.loop.bandwidth > MAX_REALISABLE_BW:
            return "loop bandwidth beyond the locking envelope"
    return ""


def evaluate(design: DesignParams, frac_bit: int | None) -> dict:
    """Analytic targets for one design point."""
    budget = noise.pd_noise_floor(design, frac_bit)
    total = noise.output_phase_noise(design, F_GRID, frac_bit)
    jitter = noise.integrated_jitter(F_GRID, total, design.f_ckv, F_LO, F_HI)

    f = design.fdvpd
    p_fdvpd = noise.fdvpd_power(design, f.i_ramp, f.comparator_noise,
                                f.c_sar)["total"]
    p_total = p_fdvpd + P_REST_W

    terms = {"ramp": budget.ramp, "ktc": budget.ktc,
             "comparator": budget.comparator, "quantisation": budget.quantisation}
    out = {
        "L_ramp_dbc": 10 * math.log10(budget.ramp) if budget.ramp > 0 else None,
        "L_ktc_dbc": 10 * math.log10(budget.ktc) if budget.ktc > 0 else None,
        "L_comparator_dbc": 10 * math.log10(budget.comparator),
        "L_quantisation_dbc": 10 * math.log10(budget.quantisation),
        "L_pd_total_dbc": 10 * math.log10(budget.pd_total),
        "dominant_pd_term": max(terms, key=terms.get),
        "jitter_s": jitter,
        "p_fdvpd_w": p_fdvpd,
        "p_total_w": p_total,
        "fom_db": noise.fom(jitter, p_total),
    }
    for probe in PROBES:
        key = f"L_out_{probe:.0f}hz_dbc".replace("000000hz", "mhz") \
            .replace("000hz", "khz")
        out[key] = 10 * math.log10(float(np.interp(probe, F_GRID, total)))
    return out


def derived(design: DesignParams) -> dict:
    f = design.fdvpd
    return {
        "f_ckv_hz": design.f_ckv,
        "f_pd_hz": design.f_pd,
        "t_pd_s": design.t_pd,
        "i_ramp_a": f.i_ramp,
        "dac_lsb_s": design.dac_lsb_time,
        "dac_lsb_v": design.dac_lsb_volt,
        "dac_fs_v": design.dac_full_scale_volt,
        "adc_fs_v": f.adc_full_scale,
        "adc_lsb_s": f.adc_lsb_time,
        "effective_bits": design.effective_bits,
    }


# --------------------------------------------------------------------------
#  Simulator cross-check
# --------------------------------------------------------------------------
def cross_check(design: DesignParams, seed: int, n_cycles: int) -> dict:
    """Run the event-driven simulator on one design and report what it did."""
    try:
        res = FdvPll(design, seed=seed).run(n_cycles)
    except Exception as exc:                      # a design the sim cannot run
        return {"sim_jitter_s": None, "sim_cycle_slips": None,
                "sim_saturated_frac": None, "sim_locked": False,
                "sim_note": type(exc).__name__}
    sat = float(res.saturated[res.discard:].mean())
    slips = res.cycle_slips
    return {
        "sim_jitter_s": res.jitter(F_LO, F_HI),
        "sim_cycle_slips": slips,
        "sim_saturated_frac": sat,
        "sim_locked": bool(sat < 1e-3 and slips == 0),
        "sim_note": "",
    }


# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=16384,
                    help="design points; rounded up to a power of two (Sobol)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--profile", choices=("wide", "realisable"),
                    default="realisable",
                    help="'realisable' also applies the two empirical lock "
                         "constraints; 'wide' keeps the full space and lets "
                         "you see where the analytic model stops being safe")
    ap.add_argument("--sim", type=int, default=0,
                    help="how many feasible rows to cross-check with the "
                         "event-driven simulator (~6 s each)")
    ap.add_argument("--sim-cycles", type=int, default=1 << 16)
    ap.add_argument("--out", type=pathlib.Path,
                    default=DATA / "design_space.csv")
    args = ap.parse_args()

    n_pow2 = max(1, math.ceil(math.log2(max(args.n, 2))))
    n = 2 ** n_pow2
    if n != args.n:
        print(f"  rounding {args.n} up to {n} to keep the Sobol balance")

    t0 = time.time()
    cols = sample_space(n_pow2, args.seed)
    print(f"  sampled {n} points over {len(DIMS)} parameters "
          f"({args.profile} profile)")

    rows: list[dict] = []
    n_feasible = 0
    for i in range(n):
        row = {k: v[i] for k, v in cols.items()}
        design = build_design(row)
        frac_bit = None if row["frac_bit"] < 0 else int(row["frac_bit"])

        record = {k: (int(v) if isinstance(v, (np.integer,)) else float(v))
                  for k, v in row.items()}
        record["integer_n"] = frac_bit is None
        record.update(derived(design))
        margins = capture_margin(design)
        record.update(margins)

        reason = check(design, margins["pd_capture_margin"], args.profile)
        record["feasible"] = not reason
        record["reject_reason"] = reason
        if reason:
            record.update({k: None for k in (
                "L_ramp_dbc", "L_ktc_dbc", "L_comparator_dbc",
                "L_quantisation_dbc", "L_pd_total_dbc", "dominant_pd_term",
                "jitter_s", "p_fdvpd_w", "p_total_w", "fom_db")})
            for probe in PROBES:
                key = f"L_out_{probe:.0f}hz_dbc".replace("000000hz", "mhz") \
                    .replace("000hz", "khz")
                record[key] = None
        else:
            n_feasible += 1
            record.update(evaluate(design, frac_bit))

        record.update({"sim_jitter_s": None, "sim_cycle_slips": None,
                       "sim_saturated_frac": None, "sim_locked": None,
                       "sim_note": ""})
        rows.append(record)

        if (i + 1) % 4096 == 0:
            print(f"    {i+1}/{n} evaluated ({time.time()-t0:.0f} s)")

    print(f"  {n_feasible}/{n} feasible ({n_feasible/n*100:.1f} %)")

    # ---- simulator cross-check on a random feasible subset ---------------
    if args.sim > 0:
        rng = np.random.default_rng(args.seed)
        candidates = [i for i, r in enumerate(rows) if r["feasible"]]
        chosen = rng.choice(candidates, size=min(args.sim, len(candidates)),
                            replace=False)
        print(f"  cross-checking {len(chosen)} rows with the simulator "
              f"({args.sim_cycles} cycles each)")
        for k, i in enumerate(chosen):
            row = rows[i]
            design = build_design({d: row[d] for d in DIMS})
            rows[i].update(cross_check(design, args.seed, args.sim_cycles))
            j_a, j_s = row["jitter_s"], rows[i]["sim_jitter_s"]
            tag = "locked" if rows[i]["sim_locked"] else "NOT locked"
            ratio = f"{j_s/j_a:6.3f}" if (j_s and j_a) else "   -  "
            print(f"    [{k+1}/{len(chosen)}] row {i}: analytic "
                  f"{j_a*1e15:8.1f} fs, sim "
                  f"{(j_s*1e15 if j_s else float('nan')):8.1f} fs  "
                  f"ratio {ratio}  {tag}")

    # ---- write ----------------------------------------------------------
    DATA.mkdir(exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {_show(args.out)}  "
          f"({len(rows)} rows x {len(fields)} columns)")

    meta = {
        "generator": "scripts/make_dataset.py",
        "fdvadpll_version": __version__,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "seed": args.seed,
        "profile": args.profile,
        "n_rows": n,
        "n_feasible": n_feasible,
        "labels_from": "analytic s-domain model (fdvadpll.noise), Sec. II-B "
                       "equations (1)-(13) plus the closed-loop composite",
        "sampling": "scrambled Sobol, log-uniform where marked",
        "fixed": {
            "f_ref_hz": F_REF,
            "jitter_band_hz": [F_LO, F_HI],
            "ramp_headroom_v": RAMP_HEADROOM_V,
            "i_ramp_bounds_a": [I_RAMP_MIN, I_RAMP_MAX],
            "p_rest_w": P_REST_W,
            "p_integer_n": P_INTEGER_N,
            "min_capture_margin": MIN_CAPTURE_MARGIN if args.profile ==
                                  "realisable" else None,
            "max_realisable_bw_hz": MAX_REALISABLE_BW if args.profile ==
                                    "realisable" else None,
        },
        "continuous_ranges": {k: {"low": v[0], "high": v[1], "scale": v[2]}
                              for k, v in CONTINUOUS.items()},
        "discrete_values": {k: list(v) for k, v in DISCRETE.items()},
        "frac_bits": list(FRAC_BITS),
        "conventions": {
            "frac_bit": "-1 means an integer-N channel; see integer_n",
            "feasible": "False rows keep their design columns and have empty "
                        "targets; reject_reason says which constraint failed",
            "sim_*": "empty except on the --sim subset",
        },
        "limitations": [
            "The oscillator's power is not modelled: dco_pn_1mhz_dbc varies "
            "freely but costs nothing, so p_total_w and fom_db understate the "
            "true cost of a quiet DCO. Treat fom_db as a detector figure of "
            "merit.",
            "p_total_w is the computed FDVPD power plus a fixed allowance for "
            "everything else, taken from the prototype's measured breakdown.",
            "Labels are analytic. The event-driven simulator agrees to a few "
            "per cent on the default design point, but a design that does not "
            "lock has no meaningful analytic jitter -- use --sim to see where.",
            "The I-DAC is not dithered, so deep fractional channels carry a "
            "DAC-quantisation floor the analytic model does not include.",
        ],
    }
    meta_path = args.out.with_suffix(".meta.json")
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"  wrote {_show(meta_path)}")
    print(f"\n  done in {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
