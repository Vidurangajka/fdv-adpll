"""Event-driven simulator for the FDV-domain fractional-N DPLL / ADPLL.

Phase bookkeeping
-----------------
Everything is referred to the *divided* clock CKVd, which is what the FDVPD
actually sees, and converted to output (CKV) cycles at the end.

    Phi_R[k]  ideal accumulated CKVd phase at reference edge k, = k * FCW_pd.
              Kept as an exact integer + fractional pair so that deep
              fractional channels (2**-16) stay bit-accurate over 10**6 cycles.
    n_edge    index of the first CKVd edge strictly after REF -- the edge the
              clock gate passes to the ramp (CKVdg).
    Phi_V[k]  actual CKVd phase at REF edge k = n_edge - dt/T_PD.

The phase error the loop acts on is ``Phi_R - Phi_V``, which splits cleanly into

    integer part    n_slip = ceil(Phi_R) - n_edge       (from the counter)
    fractional part -err_time / T_PD                    (from the FDVPD)

Writing it this way means the simulation never has to special-case the
wrap that happens when the fractional accumulator rolls over: the counter
resolves it, exactly as it does in hardware.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np

from . import dsp
from .blocks import Fdvpd
from .calib import DacInlCalibration, KdcoCalibration, PdGainCalibration
from .dco import Dco
from .dlf import DigitalLoopFilter, FrequencyLockLoop
from .params import DesignParams

__all__ = ["FdvPll", "SimResult"]


# --------------------------------------------------------------------------
@dataclass
class SimResult:
    """Output of one simulation run.  All arrays are one sample per REF cycle."""

    design: DesignParams
    tau: np.ndarray             # output edge timing error [s]
    phi_e: np.ndarray           # measured phase error [output cycles]
    code: np.ndarray            # signed SAR code
    dac_code: np.ndarray
    dt: np.ndarray              # REF -> CKVdg delay [s]
    ntw: np.ndarray             # normalised tuning word
    f_dco: np.ndarray           # instantaneous DCO frequency [Hz]
    saturated: np.ndarray
    n_slip: np.ndarray
    discard: int
    fll_disabled_at: int | None = None
    gain_history: np.ndarray | None = None
    kdco_history: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    # -- convenience -------------------------------------------------------
    @property
    def f_ref(self) -> float:
        return self.design.f_ref

    @property
    def f_ckv(self) -> float:
        return self.design.f_ckv

    @property
    def tau_settled(self) -> np.ndarray:
        return self.tau[self.discard:]

    @property
    def cycle_slips(self) -> int:
        """Number of one-period jumps in the output edge selection after lock."""
        x = self.tau_settled
        if len(x) < 2:
            return 0
        return int(np.sum(np.abs(np.diff(x)) > 0.5 * self.design.t_pd))

    def phase_noise(self, nperseg: int | None = None):
        """``(f, L)`` with ``L`` linear, sampled at ``f_REF`` (Nyquist f_REF/2)."""
        return dsp.phase_noise_from_time_error(self.tau_settled, self.f_ref,
                                               self.f_ckv, nperseg=nperseg)

    def jitter(self, f_lo: float = 10e3, f_hi: float = 40e6,
               nperseg: int | None = None) -> float:
        f, L = self.phase_noise(nperseg)
        f_hi = min(f_hi, f[-1])
        return dsp.rms_jitter(f, L, self.f_ckv, f_lo, f_hi)

    def spurs(self, n_fft: int | None = None):
        return dsp.find_spurs(self.tau_settled, self.f_ref, self.f_ckv, n_fft)

    def worst_spur(self, f_lo: float = 10e3, f_hi: float | None = None):
        """``(f, dBc)`` of the largest discrete tone in the band."""
        f, dbc = self.spurs()
        f_hi = f[-1] if f_hi is None else f_hi
        m = (f >= f_lo) & (f <= f_hi)
        if not m.any():
            return float("nan"), float("-inf")
        i = int(np.argmax(dbc[m]))
        return float(f[m][i]), float(dbc[m][i])

    def fractional_spur(self, n_harmonic: int = 1):
        """Level of the ``n``-th harmonic of the fractional frequency."""
        frac = self.design.fcw - round(self.design.fcw)
        if abs(frac) < 1e-12:
            return float("nan"), float("-inf")
        f_frac = abs(frac) * self.f_ref
        return dsp.spur_level(self.tau_settled, self.f_ref, self.f_ckv,
                              n_harmonic * f_frac)

    def power(self) -> float:
        return self.design.power.p_total_measured

    def fom(self, **kw) -> float:
        from .noise import fom
        return fom(self.jitter(**kw), self.power())

    def summary(self) -> str:
        j = self.jitter()
        fs, sd = self.worst_spur()
        lines = [
            f"  f_CKV              : {self.f_ckv/1e9:.6f} GHz "
            f"(FCW = {self.design.fcw:.6f})",
            f"  rms jitter 10k-40M : {j*1e15:.1f} fs",
            f"  FoM                : {self.fom():.1f} dB "
            f"@ {self.power()*1e3:.1f} mW",
            f"  worst spur         : {sd:.1f} dBc @ {fs/1e3:.1f} kHz",
            f"  cycle slips        : {self.cycle_slips}",
        ]
        if self.fll_disabled_at is not None:
            lines.append(f"  FLL off at         : cycle {self.fll_disabled_at}"
                         f" ({self.fll_disabled_at/self.f_ref*1e6:.2f} us)")
        frac = self.design.fcw - round(self.design.fcw)
        if abs(frac) > 1e-12:
            ff, fd = self.fractional_spur()
            lines.append(f"  fractional spur    : {fd:.1f} dBc @ "
                         f"{ff/1e3:.1f} kHz")
        return "\n".join(lines)


# --------------------------------------------------------------------------
class FdvPll:
    """Top-level DPLL / ADPLL.

    Parameters
    ----------
    design
        Design point (see :mod:`fdvadpll.params`).
    dco_mode
        ``'ideal'`` reproduces the paper; ``'banks'`` switches on the ADPLL
        tuning-bank model with sigma-delta dithering of the tracking bank.
    integer_path
        ``'fll'``   counter-based frequency-lock path that is switched off once
                    the ADC stops saturating -- the paper's Fig. 2.
        ``'phase'`` the counter's integer phase error is fed to the loop at all
                    times (the classic ADPLL variable-phase accumulator).
        ``'none'``  fractional path only; useful for isolating PD behaviour.
    calibrate
        Any of ``'gain'``, ``'inl'``, ``'kdco'``.
    f_init_offset
        Free-running centre-frequency error [Hz] the loop has to acquire, i.e.
        an untrimmed oscillator.  It shifts ``f_center`` for the whole run.
    """

    def __init__(self, design: DesignParams, seed: int = 0,
                 dco_mode: str = "ideal", integer_path: str = "fll",
                 calibrate: tuple = (), dco_noise: bool = True,
                 ref_noise: bool = True, pd_noise: bool = True,
                 iir_poles: tuple = (), f_init_offset: float = 0.0,
                 cal_mu: dict | None = None):
        if integer_path not in ("fll", "phase", "none"):
            raise ValueError("integer_path must be 'fll', 'phase' or 'none'")
        if not pd_noise:
            # all three thermal terms, not just the two with a device knob
            design = replace(design, fdvpd=replace(
                design.fdvpd, comparator_noise=0.0, gamma=0.0,
                ktc_enabled=False, flicker_corner=0.0))
        self.design = design
        self.seed = seed
        self.dco_mode = dco_mode
        self.integer_path = integer_path
        self.calibrate = set(calibrate)
        self.dco_noise = dco_noise
        self.ref_noise = ref_noise
        self.pd_noise = pd_noise
        self.iir_poles = iir_poles
        self.f_init_offset = f_init_offset
        self.cal_mu = cal_mu or {}

    # ------------------------------------------------------------------
    def run(self, n_cycles: int = 1 << 17, discard: int | None = None,
            progress: bool = False) -> SimResult:
        d = self.design
        rng = np.random.default_rng(self.seed)
        n = int(n_cycles)
        discard = n // 8 if discard is None else int(discard)

        t_ref_period = 1.0 / d.f_ref
        t_pd_ideal = 1.0 / d.f_pd
        fcw_pd = d.fcw_pd
        fcw_int = int(math.floor(fcw_pd))
        fcw_frac = fcw_pd - fcw_int
        div = d.fb_div

        # ---- noise sources ------------------------------------------------
        if self.ref_noise:
            ref_jit = dsp.time_error_from_phase_noise(
                n, d.f_ref, d.f_ref,
                lambda f: 10.0 ** (d.ref.phase_noise(f) / 10.0), rng)
        else:
            ref_jit = np.zeros(n)

        n_edges = int(n * fcw_pd) + 4096
        dco = Dco(d, rng, n_edges, mode=self.dco_mode, noise=self.dco_noise)
        # A free-running centre error the loop has to tune out, not just a
        # seed value: every tuning update recomputes f_DCO from f_center, so
        # anything written with set_frequency alone would be discarded on the
        # first cycle and the run would start on frequency after all.
        dco.f_center += self.f_init_offset
        dco.set_frequency(dco.f_center)

        pd = Fdvpd(d, rng)
        # V_OS margin, expressed as a fraction of one PD period.  Folding it
        # into the reference accumulator (rather than adding it to dt) is what
        # keeps the predicted delay inside one period for every fractional
        # word: the accumulator is simply evaluated at Psi = Phi_R - delta.
        delta = (d.fdvpd.t_offset / t_pd_ideal) % 1.0
        mean_t_on = (d.fdvpd.t_offset if abs(d.fcw - round(d.fcw)) < 1e-12
                     else 0.5 * t_pd_ideal)
        pd.prime(n, mean_t_on=mean_t_on)
        dlf = DigitalLoopFilter(d, iir_poles=self.iir_poles)
        fll = FrequencyLockLoop(d) if self.integer_path == "fll" else None

        cal_g = (PdGainCalibration(mu=self.cal_mu.get("gain", 2e-3))
                 if "gain" in self.calibrate else None)
        cal_inl = (DacInlCalibration(n_segments=2 ** d.fdvpd.dac_msb_thermo_bits,
                                     mu=self.cal_mu.get("inl", 2e-4),
                                     dac_bits=d.fdvpd.dac_bits,
                                     bin_bits=d.fdvpd.dac_bits
                                     - d.fdvpd.dac_msb_thermo_bits)
                   if "inl" in self.calibrate else None)
        cal_k = (KdcoCalibration(dco.kdco_hat, mu=self.cal_mu.get("kdco", 5e-3))
                 if "kdco" in self.calibrate else None)

        # ---- output buffers -----------------------------------------------
        tau = np.empty(n)
        phi_e_a = np.empty(n)
        code_a = np.empty(n, dtype=np.int32)
        dac_a = np.empty(n, dtype=np.int32)
        dt_a = np.empty(n)
        ntw_a = np.empty(n)
        fdco_a = np.empty(n)
        sat_a = np.zeros(n, dtype=bool)
        slip_a = np.zeros(n, dtype=np.int32)

        # ---- state ---------------------------------------------------------
        n_int_acc = 0          # integer part of Phi_R
        frac_acc = 0.0         # fractional part of Phi_R, in [0, 1)
        n_edge = 0             # index of the last generated CKVd edge
        t_edge = 0.0           # its time
        t_pd = d.t_pd
        ntw = 0.0
        periods = dco.periods
        step_pd = pd.step

        for k in range(n):
            t_ref = k * t_ref_period + ref_jit[k]

            # -- reference accumulator, evaluated at Psi = Phi_R - delta -----
            #    (exact: the integer and fractional parts are kept separate so
            #     a 2**-16 fractional word stays bit-accurate over 10**6 cycles)
            psi_frac = frac_acc - delta
            psi_int = n_int_acc
            if psi_frac < 0.0:
                psi_frac += 1.0
                psi_int -= 1
            n_exp = psi_int + (1 if psi_frac > 0.0 else 0)
            t_frac_norm = (n_exp - psi_int) - psi_frac       # in [0, 1)

            # -- advance CKVd edges up to the one the counter selects --------
            m = n_exp - n_edge
            if m > 0:
                t_new = t_edge + np.cumsum(periods(m))
                # how many of them landed at or before REF -> counter reading
                count = n_edge + int(np.searchsorted(t_new, t_ref, "right"))
                t_edge = t_new[-1]
                n_edge = n_exp
            else:
                count = n_edge
            dt = t_edge - t_ref
            # Integer part of the phase error.  The ``-1`` applies only while
            # the predicted edge really does follow REF: the phase detector
            # measures ``dt`` continuously *through* zero, so once the edge
            # slips just ahead of REF the counter would otherwise report a
            # whole extra cycle that the residue has already accounted for.
            # That double count is invisible in an integer channel -- the ramp
            # window never goes near zero there -- but in a fractional one the
            # window sweeps through zero every sawtooth period and injects a
            # full-cycle kick each time.
            n_slip = n_exp - count - (1 if dt > 0.0 else 0)

            # -- phase detection ---------------------------------------------
            tf = t_frac_norm
            if cal_g is not None:
                tf = cal_g.correct(tf)
            if cal_inl is not None:
                tf = cal_inl.correct(tf, int(round(tf * pd.dac.n_codes)))
            s = step_pd(dt, tf)

            phi_ckvd = -s.error_time / t_pd
            if self.integer_path != "none":
                phi_ckvd += n_slip
            phi_e = div * phi_ckvd                    # output CKV cycles

            # -- loop ---------------------------------------------------------
            ntw = dlf.step(phi_e)
            if fll is not None:
                ntw += fll.step(count, s.saturated, k)
            dco.set_normalised_tuning(ntw)

            # -- background calibration ---------------------------------------
            # Two gates, both necessary.
            #
            # Not while the loop is still gearing: during acquisition the phase
            # error is dominated by the frequency transient, which correlates
            # with the fractional sawtooth the gain LMS regresses on.  The
            # estimate then moves on a gradient that has nothing to do with the
            # mismatch, reaches its clamp, and strands the loop out of range
            # where no later sample can correct it.
            #
            # Never on a railed sample: a saturated ADC reports its end code,
            # not the error, so the "gradient" is whatever the clip happened to
            # be.  Both failures bite hardest in the deep fractional channels,
            # where the sawtooth is slow enough for the drift to build up.
            if not dlf.gearing:
                if not s.saturated:
                    if cal_g is not None:
                        cal_g.update(phi_e, t_frac_norm)
                    if cal_inl is not None:
                        cal_inl.update(phi_e, s.dac_code)
                elif cal_g is not None:
                    # railed: no gradient, but the leak still runs so a bad
                    # excursion decays instead of sticking
                    cal_g.relax()
            if cal_k is not None:
                dco.kdco_hat = cal_k.update(dco.otw_trk, count, d.f_ref, div,
                                            dco.coarse_hz)

            # -- record --------------------------------------------------------
            tau[k] = t_edge - n_edge * t_pd_ideal
            phi_e_a[k] = phi_e
            code_a[k] = s.code
            dac_a[k] = s.dac_code
            dt_a[k] = dt
            ntw_a[k] = ntw
            fdco_a[k] = dco.f_dco
            sat_a[k] = s.saturated
            slip_a[k] = n_slip

            # -- advance the reference accumulator ------------------------------
            n_int_acc += fcw_int
            frac_acc += fcw_frac
            if frac_acc >= 1.0:
                frac_acc -= 1.0
                n_int_acc += 1

            if progress and (k & 0xFFFF) == 0xFFFF:
                print(f"    {k+1}/{n} cycles", flush=True)

        return SimResult(
            design=d, tau=tau, phi_e=phi_e_a, code=code_a, dac_code=dac_a,
            dt=dt_a, ntw=ntw_a, f_dco=fdco_a, saturated=sat_a, n_slip=slip_a,
            discard=discard,
            fll_disabled_at=(fll.disabled_at if fll is not None else None),
            gain_history=(np.array(cal_g.history) if cal_g else None),
            kdco_history=(np.array(cal_k.history) if cal_k else None),
            meta={"dco_mode": self.dco_mode, "integer_path": self.integer_path,
                  "calibrate": sorted(self.calibrate), "seed": self.seed,
                  "dac_inl_lsb": float(np.max(np.abs(pd.dac.inl))),
                  "dac_dnl_lsb": float(np.max(np.abs(pd.dac.dnl)))},
        )
