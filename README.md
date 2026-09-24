# fdv-adpll

Behavioural model of a fractional-N digital PLL whose phase error is quantised
in the **fully differential voltage domain**, plus an all-digital (ADPLL)
extension with background calibration.

Reference implementation of

> L. Wu, T. Burger, P. Schoenle and Q. Huang, *"A Power-Efficient Fractional-N
> DPLL With Phase Error Quantized in Fully Differential-Voltage Domain,"*
> IEEE J. Solid-State Circuits, vol. 56, no. 4, pp. 1254-1264, April 2021.

The package carries two independent models of the same design:

* **`fdvadpll.noise`** -- the s-domain analysis of Sec. II-B, equations (1) to
  (13), plus the closed-loop composite.
* **`fdvadpll.pll`** -- an event-driven simulator that advances real CKVd edges
  through the real blocks (segmented I-DAC, ramp, SAR ADC, counter, loop
  filter, oscillator).

They share nothing but the design parameters, so their agreement is the main
evidence that either one is right.  On the default design point they land
within a few per cent of each other and on the paper's measured numbers:

| channel | simulated | analytic | measured |
|---|---|---|---|
| integer-N, rms jitter 10 kHz - 40 MHz | 81 fs | 84 fs | 82 fs |
| integer-N, L(110 kHz) | -116 dBc/Hz | -116 dBc/Hz | -116 dBc/Hz |
| fractional-N (FCW = 20.9688 x 2) | 94 fs | 97 fs | 101 fs |

## What this adds to the paper

The architecture, the equations and the measurements are Wu et al.'s. What is
new here is the *method* and what the method found. The full write-up is
**[Cross-Validating Two Models of a Differential-Voltage-Domain Fractional-N
DPLL](https://claude.ai/artifact/KtZU5YPiCFHPCyr2bFgeLV)**; in brief:

1. **Two models, one parameter set, as a standing check.** `noise.py` and
   `pll.py` share nothing but `params.py`, so disagreement is a signal available
   at every change rather than a single comparison at the end.
2. **Fractional-N is where the bugs are, structurally.** Six of the eight
   defects cross-validation found are invisible in integer-N. Not because that
   code was newer: an integer channel holds the internal quantities constant
   that a defect needs held. §5.1 of the report names the three mechanisms.
   Before repair the model matched the paper's *integer-N* numbers to 3 % while
   producing 51.6 ns in the fractional channel it was meant to reproduce at
   101 fs.
3. **Analytic PLL labels need two warnings attached.** The s-domain model has no
   representation of lock at all, and where the loop does lock its jitter is
   optimistic by a factor that grows with loop bandwidth. Both are quantified
   under [Design-space dataset](#design-space-dataset), along with a derived
   capture-margin feature that makes the analytic model lock-aware enough to
   filter on.
4. **The model is design entry, not documentation.** `silicon/` generates a
   sky130A specification from the same code that reproduces the measurements,
   with tolerance limits measured at loop level. See
   [`silicon/README.md`](silicon/README.md).
5. **Two silent verification failures, reported rather than quietly fixed** — a
   circuit that met every performance target while biased above its supply, and
   a tolerance sweep that quoted the first *failing* value as the permitted one.
   Both are in `silicon/README.md`; both now have a guard.

## Install

```bash
pip install -e ".[dev]"      # numpy, scipy, pytest, matplotlib
```

## Quick start

```python
from fdvadpll import FdvPll, measured_fit_design

pll = FdvPll(measured_fit_design(), seed=1)
res = pll.run(1 << 17)
print(res.summary())
```

```
  f_CKV              : 3.360000 GHz (FCW = 42.000000)
  rms jitter 10k-40M : 80.2 fs
  FoM                : -252.3 dB @ 9.2 mW
  worst spur         : -72.8 dBc @ 10.5 kHz
  cycle slips        : 0
  FLL off at         : cycle 276 (3.45 us)
```

For the analytic model alone:

```python
import numpy as np
from fdvadpll import measured_fit_design, noise

d = measured_fit_design()
f = np.logspace(3, 7.6, 2000)
L = noise.output_phase_noise(d, f, frac_bit=5)
print(noise.integrated_jitter(f, L, d.f_ckv) * 1e15, "fs")
print(noise.pd_noise_floor(d, frac_bit=5))
```

## Figures

```bash
python scripts/run_all.py            # everything, a few minutes
python scripts/run_all.py --fast     # analytic figures only, ~20 s
```

Output lands in `results/`:

| script | what it shows |
|---|---|
| `design_summary.py` | the design point and its noise budget, as text |
| `fig05_noise_vs_power.py` | best in-band FDVPD floor versus detector power |
| `fig11_phase_noise.py` | simulated vs analytic vs measured, both channels |
| `fig12_noise_breakdown.py` | which source contributes what, across the band |
| `fig15_spur_vs_fraction.py` | fractional spur versus asserted fractional bit |
| `fig17_power_and_fom.py` | the 9.2 mW breakdown and the jitter-power trade |
| `compare_architectures.py` | this detector versus a MASH-dithered divider |
| `calibration_demo.py` | acquisition and the background gain calibration |
| `make_dataset.py` | design-space dataset, see below |

## Design-space dataset

```bash
python scripts/make_dataset.py                        # 16384 points, ~30 s
python scripts/make_dataset.py --n 65536 --sim 64     # bigger, + sim check
python scripts/make_dataset.py --profile wide         # no lock constraints
```

Sobol-samples 18 design parameters, evaluates the analytic model at each point,
and writes `data/design_space.csv` plus a `.meta.json` recording the ranges,
the fixed choices and the limitations. 53 columns: design parameters and
derived quantities as features; per-term in-band noise, the output profile at
four probe offsets, jitter, power and FoM as targets.

```python
import csv
rows = [r for r in csv.DictReader(open("data/design_space.csv"))
        if r["feasible"] == "True"]
```

Generation is deterministic given `--seed`, so the 12 MB CSV does not need to
be kept — regenerate it in about 30 s. The shipped default (16384 points,
`realisable` profile, seed 0):

| | |
|---|---|
| feasible rows | 7062 / 16384 (43 %) |
| analytic jitter, feasible rows | 45 / 128 / 394 fs (p5 / median / p95) |
| FoM range | −265 to −230 dB |
| dominant in-band term | ramp current 46 %, comparator 32 %, kT/C 20 %, quantisation 2 % |

### Read this before training anything on it

**The analytic model has no notion of lock.** `noise.py` will happily return a
jitter figure for a design whose detector would be railed from the first cycle,
because nothing in the s-domain analysis knows the ADC has ends. The
`realisable` profile therefore adds two constraints that do *not* follow from
the model — a minimum capture margin and a bandwidth ceiling — both measured by
cross-checking Sobol samples against the simulator. `--profile wide` drops them
if you want to learn where the envelope is.

**The labels are optimistic, and how much depends on loop bandwidth.**
Cross-checking 48 feasible rows with the event-driven simulator: 38 locked
(79 %), and on those the simulated jitter ran

| loop bandwidth | median sim / analytic |
|---|---|
| < 200 kHz | 1.50× |
| 200 kHz – 500 kHz | 4.65× |
| ≥ 500 kHz | 4.54× |

12 of 38 within 1.5×, 25 of 38 within 3×. The likely cause is that
`loop_transfer` is a continuous-time approximation of a loop that is actually
sampled at `f_REF`; the prototype sits at 500 kHz and agrees to 3 %, but that
agreement does not survive moving the other parameters far off nominal. It is
*not* driven by the fractional word — deep fractional channels were, if
anything, better. Treat the analytic labels as a fast, smooth, systematically
optimistic surrogate rather than ground truth, and use `--sim` to re-measure
the bias on whatever region you care about. (38 points is a small sample;
these are indicative.)

**Other caveats** are listed in the `.meta.json`, notably that the oscillator's
power is not modelled, so `fom_db` is a detector figure of merit rather than a
whole-synthesiser one.

## Design points

`params.py` offers three:

* **`default_design()`** -- the published equations and published device values,
  nothing else.  It comes out 4 to 7 dB optimistic in band, which Sec. V says
  it should: the s-domain analysis omits ramp-generator flicker noise.
* **`measured_fit_design()`** -- adds exactly one parameter beyond the paper,
  `fdvpd.flicker_corner = 520 kHz`, the term Sec. V itself identifies as
  missing.  A least-squares fit was free to add a white excess term too and
  converged to 0.0 dB, so the published thermal equations need no fudge.
* **`fractional_design(bit)`** -- a near-integer channel with only fractional
  bit `bit` asserted.  `bit=5` is the Fig. 11 channel; sweeping 2 to 16 gives
  Fig. 15.

## What the model does and does not cover

**Covered.**  Segmented I-DAC INL/DNL and code-dependent settling; ramp gain
error and 2nd/3rd-order ramp non-linearity; kT/C sampling noise; comparator
noise; SAR CDAC mismatch; counter-based frequency acquisition and its handover
to the phase detector; type-II loop filtering with gear shifting and optional
extra IIR poles; oscillator period jitter with a 1/f^3 corner; the ADPLL tuning
banks with sigma-delta dither; background gain, INL and K_DCO calibration.

**Known limits, all of them real rather than bookkeeping:**

* **The I-DAC is not dithered.**  `T_frac` is rounded onto the 10 b grid, so at
  fractional bit 11 and beyond the word advances by half an LSB per cycle and
  the rounding error becomes a slow in-band triangle rather than white noise.
  The 581 fs DAC LSB therefore sets a floor on fractional phase accuracy in the
  deepest channels.  Hardware answers this with dither; the model does not.
* **The gain calibration assumes its regressor survives the loop.**  It
  correlates the measured error against the fractional sawtooth.  Once that
  sawtooth falls inside the loop bandwidth, what reaches the detector is a
  phase-rotated version, and past 90 degrees of rotation the gradient changes
  sign.  The clamp and the leak keep that from costing lock, but in the deepest
  channels the loop can be inert rather than converging.  `calib.py` explains
  why running the calibration at a reduced loop bandwidth is the fix.
* **K_DCO calibration needs excitation.**  In lock the tuning word barely
  moves, so essentially all of the useful updates happen during acquisition.
  The counter is also a coarse frequency meter: one window resolves
  `f_REF*div/window`, so the window has to be tens of thousands of reference
  cycles before a 12 kHz/LSB bank is observable at all.
* **A positive DAC gain error costs lock in some deep fractional channels**,
  from a few tenths of a percent upwards: the DAC asks for more than one PD
  period at the top of the sawtooth, so the top of the range is unreachable and
  the detector rails once per sawtooth period.  Which channels break is
  deterministic (identical across seeds) but not monotonic in the fractional
  bit -- it depends on where the sawtooth harmonics land relative to the loop
  resonance.  The effect is physical, not a modelling artefact, and the
  calibration cannot help once the loop is already out of lock.

## Tests

```bash
pytest                       # 266 tests, ~85 s
pytest -m "not slow"         # skip the 20 full-length simulation runs, ~23 s
```

The suite is structured around the paper rather than around the code: derived
design numbers are checked against the values the text quotes, each equation
against its closed form, each block against the imperfection it is there to
model, and the simulator against the analytic model.

## Layout

```
fdvadpll/
  params.py   design points and every derived number, with the paper's values
  noise.py    s-domain model, equations (1)-(13) and the closed-loop composite
  blocks.py   I-DAC, ramp generator, SAR ADC, assembled phase detector
  dco.py      oscillator: ideal tuning law and the ADPLL tuning banks
  dlf.py      type-II digital loop filter and the frequency-lock path
  sdm.py      MASH modulators -- the dither, and the architecture compared against
  calib.py    background gain / INL / K_DCO calibration
  dsp.py      noise synthesis and spectral estimation
  pll.py      the event-driven simulator
scripts/      figure reproduction and dataset generation
results/      generated figures
data/         generated dataset + its provenance file
tests/        pytest suite
```
