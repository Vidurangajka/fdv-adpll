"""Analytic (s-domain) noise model of the FDV-domain DPLL.

Implements Sec. II-B of Wu et al., JSSC 2021, equations (1)-(13), plus the
closed-loop composite used for Fig. 5 and Fig. 12.

A note on the feedback divide-by-2
----------------------------------
The paper deliberately ignores the div-2 in the analysis (Sec. II-B).  It turns
out that only *one* of the mechanisms actually cares about it:

  * ``K_R`` (6) and ``K_SR`` (10) are Volts (or Amps) per radian of **output**
    phase.  A timing error is a timing error, so both are independent of the
    feedback divider.
  * (11) kT/C and (12) comparator noise therefore keep their published form.
  * (1) the ramp gating window scales with the PD period, so the div-2 doubles
    the duty cycle and costs +3 dB.
  * (13) quantisation is written here in terms of the *time* LSB instead of
    ``T_CKV/2**N_B``, which removes the ambiguity: 12 bit over one 595 ps PD
    period is 11 bit over one 298 ps CKV period.

Set ``include_div2=False`` to get the paper's published expressions verbatim.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .params import KT, DesignParams

__all__ = [
    "gating_duty_cycle",
    "gating_harmonics",
    "ramp_current_psd",
    "k_ramp",
    "k_slope",
    "L_ramp_current",
    "L_charge_pump_reference",
    "L_ktc",
    "L_comparator",
    "L_quantisation",
    "pd_noise_floor",
    "NoiseBudget",
    "loop_transfer",
    "output_phase_noise",
    "integrated_jitter",
    "fom",
    "fdvpd_power",
    "noise_vs_power",
]


# --------------------------------------------------------------------------
#  Ramp-generator noise -- equations (1) to (8)
# --------------------------------------------------------------------------
def gating_duty_cycle(design: DesignParams, frac_bit: int | None = None,
                      t_on_integer: float = 100e-12,
                      include_div2: bool = True) -> float:
    """Average ramp on-time as a fraction of ``T_REF``  --  equation (1).

    In a fractional-N channel with only fractional bit ``N_q`` asserted the
    window width steps through ``i/2**Nq * T`` for ``i = 0 .. 2**Nq-1`` and
    repeats every ``2**Nq`` reference cycles, giving a mean on-time

        g_bar = (1/2 - 1/2**(Nq+1)) * T

    In an integer-N channel the window is nominally zero; the paper models the
    residual with a fixed 100 ps on-time (Fig. 12 discussion), which is what
    ``frac_bit=None`` returns here.
    """
    t_period = design.t_pd if include_div2 else design.t_ckv
    t_ref = 1.0 / design.f_ref
    if frac_bit is None:
        return t_on_integer / t_ref
    g_bar = (0.5 - 0.5 ** (frac_bit + 1)) * t_period
    return g_bar / t_ref


def gating_harmonics(frac_bit: int, n_harmonics: int | None = None) -> float:
    """``sum_n |G_n|**2`` of the gating window  --  equations (3), (4), (5).

    Evaluated exactly rather than by truncating the series: for a pulse train
    of unit height the Fourier-series Parseval identity gives

        sum_n |G_n|^2 = (1/T_p) * integral g(t)^2 dt = duty cycle

    ``n_harmonics`` is accepted for the explicit finite-sum cross-check used in
    the unit tests.
    """
    m = 2 ** frac_bit
    duty = np.arange(m) / m          # on-time of each cycle, in units of T_pd
    if n_harmonics is None:
        # Parseval: mean of g(t)^2 over the repetition period.
        return float(duty.mean())
    # Explicit truncated Fourier sum, for validation only.
    total = 0.0
    for n in range(-n_harmonics, n_harmonics + 1):
        # G_n of a pulse train with per-cycle widths `duty` (units of T_REF
        # are folded out; only the shape matters for the identity check).
        k = np.arange(m)
        if n == 0:
            gn = duty.mean()
        else:
            # G_n = (1/T_p) int g(t) exp(-j2pi n t/T_p) dt with T_p = m*T_pd.
            # Substituting t = (k+u)*T_pd turns the 1/T_p and the du -> dt
            # Jacobian into 1/m and T_pd, which cancel; what is left is the
            # per-cycle integral 1/(j2pi n) * (1 - exp(-j2pi n duty_k/m)).
            phase = np.exp(-2j * np.pi * n * k / m)
            gn = np.sum(phase * (1 - np.exp(-2j * np.pi * n * duty / m))
                        / (2j * np.pi * n))
        total += abs(gn) ** 2
    return float(total)


def ramp_current_psd(design: DesignParams) -> float:
    """White current-noise PSD of the two ramp sources  --  equation (2).

    ``S_IR,n = 4kT*gamma * 2*gm``  [A^2/Hz]
    """
    f = design.fdvpd
    return 4.0 * KT * f.gamma * 2.0 * f.gm_ramp


def k_ramp(design: DesignParams) -> float:
    """Ramp feedback gain ``K_R = I_R / (pi*N)``  [A/rad]  --  equation (6)."""
    return design.fdvpd.i_ramp / (math.pi * design.n_mult)


def k_slope(design: DesignParams) -> float:
    """``K_SR = 2*I_R / (2*pi*f_CKV*C_SAR)``  [V/rad]  --  equation (10)."""
    f = design.fdvpd
    return 2.0 * f.i_ramp / (2.0 * math.pi * design.f_ckv * f.c_sar)


def L_ramp_current(design: DesignParams, frac_bit: int | None = None,
                   include_div2: bool = True) -> float:
    """In-band output phase noise from ``I_R``  --  equations (5) + (7).

    Returns a linear single-sideband density [1/Hz].
    """
    duty = gating_duty_cycle(design, frac_bit, include_div2=include_div2)
    s_gated = ramp_current_psd(design) * duty
    return 0.5 * s_gated / k_ramp(design) ** 2


def L_charge_pump_reference(design: DesignParams, i_cp: float | None = None,
                            t_on: float = 0.75e-9) -> float:
    """Charge-pump IB noise of an equivalent CPPLL  --  equation (8).

    Used only for the >10 dB comparison quoted in Sec. II-B1.  The published
    equation prints ``I_R`` in the denominator; that is a typo for the
    charge-pump current, which is what is used here.
    """
    f = design.fdvpd
    i_cp = f.i_ramp if i_cp is None else i_cp
    gm = f.gm_id * i_cp
    return (16.0 * math.pi ** 2 * design.n_mult ** 2 * KT * f.gamma * gm
            / i_cp ** 2 * t_on * design.f_ref)


# --------------------------------------------------------------------------
#  Sampling, comparator and quantisation noise -- equations (9) to (13)
# --------------------------------------------------------------------------
def L_ktc(design: DesignParams) -> float:
    """kT/C sampling noise  --  equations (9) + (11)."""
    f = design.fdvpd
    if not f.ktc_enabled:
        return 0.0
    return (2.0 * math.pi ** 2 * design.f_ckv ** 2 / design.f_ref
            * KT * f.c_sar / f.i_ramp ** 2)


def L_comparator(design: DesignParams) -> float:
    """ADC/comparator input-referred thermal noise  --  equation (12)."""
    f = design.fdvpd
    return (math.pi ** 2 * design.f_ckv ** 2 / design.f_ref
            * f.comparator_noise ** 2 * f.c_sar ** 2 / f.i_ramp ** 2)


def L_quantisation(design: DesignParams, lsb_time: float | None = None) -> float:
    """PD quantisation noise  --  equation (13), written with the time LSB.

    The published form ``(T_CKV / 2**N_B)**2`` is the square of the PD time
    LSB; using the LSB directly keeps the div-2 bookkeeping unambiguous.
    """
    lsb = design.fdvpd.adc_lsb_time if lsb_time is None else lsb_time
    return (math.pi ** 2 / 3.0 * design.f_ckv ** 2 / design.f_ref * lsb ** 2)


@dataclass
class NoiseBudget:
    """In-band output phase-noise contributions, linear [1/Hz]."""

    ramp: float
    ktc: float
    comparator: float
    quantisation: float

    @property
    def pd_total(self) -> float:
        return self.ramp + self.ktc + self.comparator + self.quantisation

    def as_dbc(self) -> dict:
        out = {}
        for k in ("ramp", "ktc", "comparator", "quantisation"):
            v = getattr(self, k)
            out[k] = 10.0 * math.log10(v) if v > 0 else float("-inf")
        out["pd_total"] = 10.0 * math.log10(self.pd_total)
        return out

    def __str__(self) -> str:
        d = self.as_dbc()
        names = {"ramp": "ramp current I_R", "ktc": "kT/C sampling",
                 "comparator": "comparator v_TN", "quantisation": "quantisation",
                 "pd_total": "FDVPD total"}
        lines = []
        for k in ("ramp", "ktc", "comparator", "quantisation", "pd_total"):
            share = getattr(self, k, self.pd_total) / self.pd_total * 100
            tag = "" if k == "pd_total" else f"   ({share:5.1f} %)"
            lines.append(f"  {names[k]:<18}: {d[k]:8.2f} dBc/Hz{tag}")
        return "\n".join(lines)


def pd_noise_floor(design: DesignParams, frac_bit: int | None = None,
                   include_div2: bool = True) -> NoiseBudget:
    """Full in-band phase-detector noise budget at the PLL output.

    ``fdvpd.excess_noise_db`` scales the three thermal terms; quantisation is
    left alone because it is set by the ADC LSB, not by device noise.
    """
    k = 10.0 ** (design.fdvpd.excess_noise_db / 10.0)
    return NoiseBudget(
        ramp=L_ramp_current(design, frac_bit, include_div2) * k,
        ktc=L_ktc(design) * k,
        comparator=L_comparator(design) * k,
        quantisation=L_quantisation(design),
    )


def pd_noise_profile(design: DesignParams, f, frac_bit: int | None = None,
                     include_div2: bool = True):
    """PD noise versus offset frequency, i.e. the white floor plus its 1/f rise.

    Returns ``(dict_of_contributions, total)``, all linear.  The flicker
    shaping ``(1 + f_corner/f)`` multiplies the thermal terms only.
    """
    f = np.asarray(f, dtype=float)
    b = pd_noise_floor(design, frac_bit, include_div2)
    fc = design.fdvpd.flicker_corner
    shape = 1.0 + fc / f if fc > 0.0 else np.ones_like(f)
    parts = {
        "ramp current": b.ramp * shape,
        "kT/C": b.ktc * shape,
        "comparator": b.comparator * shape,
        "quantisation": b.quantisation * np.ones_like(f),
    }
    return parts, sum(parts.values())


# --------------------------------------------------------------------------
#  Closed-loop transfer functions and the composite profile (Fig. 12)
# --------------------------------------------------------------------------
def loop_transfer(design: DesignParams, f, bandwidth=None, damping=None):
    """Return ``(H_lowpass, H_highpass)`` evaluated at offset frequencies ``f``.

    Continuous-time equivalent of the normalised type-II digital loop:
        omega_n = sqrt(rho)*f_REF,  zeta = alpha/(2*sqrt(rho))
    """
    alpha, rho = design.loop.alpha_rho(design.f_ref, bandwidth, damping)
    wn = math.sqrt(rho) * design.f_ref
    zeta = alpha / (2.0 * math.sqrt(rho))
    s = 2j * np.pi * np.asarray(f, dtype=float)
    den = s ** 2 + 2.0 * zeta * wn * s + wn ** 2
    h_lp = (2.0 * zeta * wn * s + wn ** 2) / den
    h_hp = s ** 2 / den
    return h_lp, h_hp


def dco_free_running_pn(design: DesignParams, f):
    """Free-running DCO phase noise L(f) [linear, 1/Hz]."""
    d = design.dco
    f = np.asarray(f, dtype=float)
    lin_1mhz = 10.0 ** (d.pn_1mhz / 10.0)
    white = lin_1mhz * (1e6 / f) ** 2
    flicker = white * (d.f_corner_flicker / f)
    return white + flicker


def output_phase_noise(design: DesignParams, f, frac_bit: int | None = None,
                       include_div2: bool = True, breakdown: bool = False):
    """Composite output phase noise L(f) [linear] -- reproduces Fig. 12.

    Returns the total, or ``(total, dict_of_contributions)`` if ``breakdown``.
    """
    f = np.asarray(f, dtype=float)
    h_lp, h_hp = loop_transfer(design, f)
    lp2 = np.abs(h_lp) ** 2
    hp2 = np.abs(h_hp) ** 2

    pd_parts, _ = pd_noise_profile(design, f, frac_bit, include_div2)
    ref_lin = 10.0 ** (design.ref.phase_noise(f) / 10.0) * design.n_mult ** 2

    parts = {"reference": ref_lin * lp2}
    parts.update({k: v * lp2 for k, v in pd_parts.items()})
    parts["DCO"] = dco_free_running_pn(design, f) * hp2
    total = sum(parts.values())
    if breakdown:
        parts["DCO (free running)"] = dco_free_running_pn(design, f)
        return total, parts
    return total


# --------------------------------------------------------------------------
#  Jitter and figure of merit
# --------------------------------------------------------------------------
def integrated_jitter(f, L_lin, f_ckv: float, f_lo: float = 10e3,
                      f_hi: float = 40e6) -> float:
    """RMS jitter [s] from a single-sideband profile, integrated ``f_lo..f_hi``.

    ``sigma_t = sqrt(2 * int L df) / (2*pi*f_ckv)``
    """
    f = np.asarray(f, dtype=float)
    L_lin = np.asarray(L_lin, dtype=float)
    m = (f >= f_lo) & (f <= f_hi)
    if m.sum() < 2:
        raise ValueError("integration band contains fewer than two points")
    ipn = 2.0 * np.trapezoid(L_lin[m], f[m])
    return math.sqrt(ipn) / (2.0 * math.pi * f_ckv)


def fom(jitter_s: float, power_w: float) -> float:
    """Jitter-power figure of merit [dB]  --  Table I footnote.

    ``FoM = 20*log10(sigma_rms/1s) + 10*log10(P/1mW)``
    """
    return 20.0 * math.log10(jitter_s) + 10.0 * math.log10(power_w / 1e-3)


def normalise_spur(spur_dbc: float, f_carrier: float,
                   f_ref_carrier: float = 3.24e9) -> float:
    """Table I's spur normalisation to a common carrier."""
    return spur_dbc + 20.0 * math.log10(f_ref_carrier / f_carrier)


# --------------------------------------------------------------------------
#  Fig. 5 : in-band noise floor versus FDVPD power
# --------------------------------------------------------------------------
def fdvpd_power(design: DesignParams, i_ramp: float, comparator_noise: float,
                c_sar: float) -> dict:
    """Power model of the FDVPD used for Fig. 5.

    ``P = P_IDAC,static + P_CSAR,dynamic + P_comparator + P_logic``

    * the static I-DAC carries the ramp branch inside it (Sec. III-C), so its
      current is modelled as ``dac_current_ratio * I_R``;
    * the differential ramp discharges ``C_SAR`` over the full-scale swing once
      per reference cycle;
    * the dynamic comparator obeys the ``30 nJ*uV^2`` energy-noise product of
      Sec. II-B5, per conversion.
    """
    p = design.power
    v_fs = design.t_pd * design.fdvpd.slew_rate
    p_dac = p.vdd * p.dac_current_ratio * i_ramp
    p_dyn = 2.0 * c_sar * v_fs * p.vdd * design.f_ref
    p_comp = (p.comparator_energy_noise / comparator_noise ** 2) * design.f_ref
    p_logic = 0.05e-3
    return {"idac": p_dac, "cdyn": p_dyn, "comparator": p_comp,
            "logic": p_logic,
            "total": p_dac + p_dyn + p_comp + p_logic}


def noise_vs_power(design: DesignParams, n_bits: int, powers,
                   frac_bit: int | None = 9, include_div2: bool = False,
                   n_grid: int = 400):
    """Best achievable FDVPD in-band noise floor for each power budget.

    For every point the slew rate is held at its nominal value (Fig. 5 fixes
    ``dv/dt``), so ``C_SAR = 2*I_R/SR``.  The split between ramp current and
    comparator power is then *optimised* to minimise the total floor, which is
    what produces the family of curves in Fig. 5.

    Returns ``(L_dbc, info)`` where ``info`` holds the optimum ``I_R``,
    ``C_SAR`` and ``v_TN`` at each power.
    """
    from dataclasses import replace as _replace

    powers = np.atleast_1d(np.asarray(powers, dtype=float))
    sr = design.fdvpd.slew_rate
    lsb_t = design.t_pd / 2 ** n_bits
    l_qn = L_quantisation(design, lsb_t)

    out = np.empty_like(powers)
    info = {"i_ramp": np.empty_like(powers), "c_sar": np.empty_like(powers),
            "v_tn": np.empty_like(powers), "feasible": np.ones_like(powers, bool)}

    pp = design.power
    for idx, budget in enumerate(powers):
        # fraction of the budget spent on the comparator
        fracs = np.linspace(0.02, 0.90, n_grid)
        best, best_state = np.inf, None
        for fr in fracs:
            p_comp = fr * budget
            v_tn = math.sqrt(pp.comparator_energy_noise * design.f_ref / p_comp)
            # remaining budget -> I_R, accounting for the dynamic term
            v_fs = design.t_pd * sr
            # P_rest = vdd*k*I_R + 2*(2 I_R/SR)*v_fs*vdd*f_ref + P_logic
            p_rest = (1.0 - fr) * budget - 0.05e-3
            if p_rest <= 0:
                continue
            coeff = (pp.vdd * pp.dac_current_ratio
                     + 4.0 * v_fs * pp.vdd * design.f_ref / sr)
            i_r = p_rest / coeff
            c_sar = 2.0 * i_r / sr
            d = _replace(design, fdvpd=_replace(design.fdvpd, c_sar=c_sar,
                                                comparator_noise=v_tn))
            tot = (L_ramp_current(d, frac_bit, include_div2)
                   + L_ktc(d) + L_comparator(d) + l_qn)
            if tot < best:
                best, best_state = tot, (i_r, c_sar, v_tn)
        if best_state is None:
            out[idx] = np.nan
            info["feasible"][idx] = False
            continue
        out[idx] = 10.0 * math.log10(best)
        info["i_ramp"][idx], info["c_sar"][idx], info["v_tn"][idx] = best_state
    return out, info
