"""Design parameters for the FDV-domain fractional-N DPLL / ADPLL.

Defaults reproduce the prototype of

    L. Wu, T. Burger, P. Schoenle and Q. Huang, "A Power-Efficient Fractional-N
    DPLL With Phase Error Quantized in Fully Differential-Voltage Domain,"
    IEEE JSSC, vol. 56, no. 4, pp. 1254-1264, April 2021.

Derived numbers (all checked against the paper text):

    f_REF     = 80 MHz,  FCW = 42  ->  f_CKV = 3.36 GHz
    fb_div    = 2        ->  the FDVPD runs at f_PD = 1.68 GHz, T_PD = 595.24 ps
                             ("FCW = 21 x 2" in the paper's notation)
    SR        = 2*I_R/C_SAR = 0.8 mV/ps        (Sec. III-A)
    C_SAR     = 550 fF                         (Sec. II-B2)
      => I_R  = SR*C_SAR/2 = 220 uA
    DAC full scale (differential) = SR*T_PD = 476 mV, i.e. 238 mV single ended,
      which sits inside the 300 mV single-ended ramp range quoted in Sec. III-B.
    DAC LSB   = T_PD/2^10 = 581 fs = 465 uV
    ADC LSB   = 118 uV = 147.5 fs              (Sec. III-A)
      => effective sub-ranging resolution = T_PD/147.5fs = 4035 ~ 12.0 bit,
         which is the "12-bit out of one CKV period" of Sec. II-B4.
    ADC range = 2^7 * 118 uV = 15.1 mV         ("slightly more than 10 mV")

The oscillator / reference phase-noise numbers are *fits* chosen to reproduce
Fig. 11 and Fig. 12; the paper does not publish the raw device data.  They are
flagged as such below.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math

# ---------------------------------------------------------------- constants --
K_BOLTZMANN = 1.380649e-23
T_ROOM = 300.0
KT = K_BOLTZMANN * T_ROOM


@dataclass
class FdvpdParams:
    """Fully-differential-voltage phase detector (Fig. 6)."""

    # --- ramp generator -------------------------------------------------
    c_sar: float = 550e-15          # F, sampling / CDAC capacitance
    slew_rate: float = 0.8e-3 / 1e-12   # V/s differential, = 0.8 mV/ps
    ramp_mismatch: float = 0.0      # fractional mismatch between the two I_R
    ramp_nl2: float = 0.0           # 2nd-order ramp coefficient  [1/V]
    ramp_nl3: float = 0.0           # 3rd-order ramp coefficient  [1/V^2]
    v_cm: float = 0.55              # V, DAC/ramp common mode
    # V_OS/SR: fixed margin keeping REF ahead of the expected CKVd edge.  It is
    # also the whole ramp window in an integer-N channel, which is why Fig. 12
    # models the integer-N I_R contribution with a 100 ps on-time.
    t_offset: float = 100e-12

    # --- 10 b segmented current DAC -------------------------------------
    dac_bits: int = 10
    dac_msb_thermo_bits: int = 4    # 4 b thermometer MSB + 6 b binary LSB
    dac_sigma_lsb: float = 0.0      # unit-element sigma, in DAC LSB
    dac_gain_err: float = 0.0       # static DAC/ramp gain mismatch (fractional)
    dac_settle_tau: float = 0.0     # s, R_D*C settling time constant

    # --- 7 b self-timed SAR ADC -----------------------------------------
    adc_bits: int = 7
    adc_lsb: float = 118e-6         # V
    adc_sigma_cap: float = 0.0      # CDAC unit-cap sigma (fractional)
    comparator_noise: float = 100e-6  # V rms, input referred (Sec. II-B3)

    # --- device noise ----------------------------------------------------
    gamma: float = 1.0              # excess noise factor
    # kT/C sampling noise, equation (9).  Unlike the ramp and comparator terms
    # it has no device parameter that can switch it off (it follows from C_SAR
    # alone), so it gets an explicit flag -- otherwise a "noiseless" phase
    # detector still carries ~123 uV rms, about one ADC LSB.
    ktc_enabled: bool = True
    gm_id: float = 8.0              # ramp-current transistor gm/ID  (Fig. 5)
    # 1/f corner of the phase-detector noise, dominated by the ramp current
    # sources.  Sec. V states this is the reason the s-domain analysis of
    # Sec. II-B under-predicts the measurement between 1 and 100 kHz; it is
    # zero in the published equations, hence zero by default here.
    flicker_corner: float = 0.0
    # Lumps PD noise that equations (1)-(13) do not cover: div-2 and clock-gate
    # jitter, REF slicer noise, supply coupling.  Zero for the pure-physics
    # model; fitted in :func:`measured_fit_design`.
    excess_noise_db: float = 0.0

    @property
    def i_ramp(self) -> float:
        """Single-ended ramp current I_R [A].  SR = 2*I_R/C_SAR."""
        return self.slew_rate * self.c_sar / 2.0

    @property
    def gm_ramp(self) -> float:
        return self.gm_id * self.i_ramp

    @property
    def adc_full_scale(self) -> float:
        """Differential input range of the SAR ADC [V]."""
        return (2 ** self.adc_bits) * self.adc_lsb

    @property
    def adc_lsb_time(self) -> float:
        """One ADC LSB expressed as a time error [s]."""
        return self.adc_lsb / self.slew_rate


@dataclass
class DcoParams:
    """Harmonic-shaped LC DCO (Fig. 9).

    ``pn_1mhz`` / ``f_corner_flicker`` are fitted to Fig. 12, they are not
    measured device data.
    """

    f_center: float = 3.36e9        # Hz, free-running centre
    kdco: float = 12e3              # Hz/LSB of the tracking bank
    kdco_err: float = 0.0           # fractional error of the *assumed* K_DCO
    pn_1mhz: float = -130.4         # dBc/Hz at 1 MHz offset (fitted to Fig. 11)
    f_corner_flicker: float = 80e3  # Hz, 1/f^3 corner (Sec. IV)
    # --- tuning banks (ADPLL extension) ---------------------------------
    pvt_bits: int = 8
    acq_bits: int = 8
    trk_bits: int = 8
    kdco_pvt: float = 4e6           # Hz/LSB
    kdco_acq: float = 300e3         # Hz/LSB
    trk_dither_bits: int = 5        # sigma-delta dithering of the tracking bank

    def period_jitter_white(self) -> float:
        """sigma of the uncorrelated period deviation that yields ``pn_1mhz``.

        L(df) = f0^3 * sigma_dT^2 / df^2   (white period jitter, 1/f^2 region)
        """
        lin = 10.0 ** (self.pn_1mhz / 10.0)
        return math.sqrt(lin * (1e6 ** 2) / self.f_center ** 3)


@dataclass
class RefParams:
    """80 MHz crystal reference.  Profile fitted to the green trace of Fig. 12."""

    f_ref: float = 80e6
    pn_floor: float = -160.0        # dBc/Hz white floor (fitted to Fig. 11)
    f_corner: float = 1e3           # Hz, 1/f^3 corner

    def phase_noise(self, f):
        """Single-sideband phase noise L(f) [dBc/Hz] of the reference."""
        import numpy as np

        f = np.asarray(f, dtype=float)
        lin = 10.0 ** (self.pn_floor / 10.0) * (1.0 + (self.f_corner / f) ** 3)
        return 10.0 * np.log10(lin)


@dataclass
class LoopParams:
    """Type-II digital loop filter, normalised ADPLL formulation.

    With ``NTW = alpha*phi_e + rho*sum(phi_e)`` and phase measured in output
    cycles, the continuous-time equivalents are

        omega_u = alpha * f_REF            (unity gain frequency, rad/s)
        omega_n = sqrt(rho) * f_REF
        zeta    = alpha / (2*sqrt(rho))
    """

    # 500 kHz is what the Fig. 11 profile implies for the prototype; it is a
    # design choice rather than a device parameter, so change it freely.
    bandwidth: float = 500e3        # Hz, target closed-loop unity-gain BW
    damping: float = 1.0
    # gear shifting: (bandwidth, damping) applied for the first N ref cycles
    gear_bandwidth: float = 6e6
    gear_damping: float = 1.0
    gear_cycles: int = 3000

    def alpha_rho(self, f_ref: float, bandwidth=None, damping=None):
        bw = self.bandwidth if bandwidth is None else bandwidth
        z = self.damping if damping is None else damping
        alpha = 2.0 * math.pi * bw / f_ref
        rho = (alpha / (2.0 * z)) ** 2
        return alpha, rho


@dataclass
class PowerParams:
    """Power model used for Fig. 5 and for the FoM bookkeeping (Fig. 17)."""

    vdd: float = 1.2                # V, everything except the DCO
    vdd_dco: float = 1.5            # V
    comparator_energy_noise: float = 30e-9 * 1e-12  # J*V^2  (30 nJ*uV^2)
    sar_cycles: int = 7
    dac_current_ratio: float = 2.4  # I_DAC,static / I_R
    p_total_measured: float = 9.2e-3
    # Fig. 17 breakdown
    share_dco_div2: float = 0.72
    share_dac_ramp: float = 0.11
    share_adc: float = 0.06
    share_dlf: float = 0.08
    share_rest: float = 0.03


@dataclass
class DesignParams:
    """Top-level design point."""

    f_ref: float = 80e6
    fcw: float = 42.0               # total FCW: f_CKV = fcw * f_ref
    fb_div: int = 2                 # divider in the feedback path (Fig. 2)

    fdvpd: FdvpdParams = field(default_factory=FdvpdParams)
    dco: DcoParams = field(default_factory=DcoParams)
    ref: RefParams = field(default_factory=RefParams)
    loop: LoopParams = field(default_factory=LoopParams)
    power: PowerParams = field(default_factory=PowerParams)

    # ---- derived ---------------------------------------------------------
    @property
    def f_ckv(self) -> float:
        return self.fcw * self.f_ref

    @property
    def t_ckv(self) -> float:
        return 1.0 / self.f_ckv

    @property
    def f_pd(self) -> float:
        """Frequency seen by the FDVPD (after the feedback divider)."""
        return self.f_ckv / self.fb_div

    @property
    def t_pd(self) -> float:
        return 1.0 / self.f_pd

    @property
    def fcw_pd(self) -> float:
        """FCW referred to the phase detector, i.e. the paper's '21 x 2' -> 21."""
        return self.fcw / self.fb_div

    @property
    def n_mult(self) -> float:
        return self.fcw

    @property
    def dac_lsb_time(self) -> float:
        return self.t_pd / (2 ** self.fdvpd.dac_bits)

    @property
    def dac_lsb_volt(self) -> float:
        return self.dac_lsb_time * self.fdvpd.slew_rate

    @property
    def dac_full_scale_volt(self) -> float:
        return self.t_pd * self.fdvpd.slew_rate

    @property
    def effective_bits(self) -> float:
        """Equivalent PD resolution over one PD period, from sub-ranging."""
        return math.log2(self.t_pd / self.fdvpd.adc_lsb_time)

    def summary(self) -> str:
        f = self.fdvpd
        rows = [
            ("f_REF", f"{self.f_ref/1e6:.1f} MHz"),
            ("FCW (total)", f"{self.fcw:.6f}"),
            ("f_CKV", f"{self.f_ckv/1e9:.5f} GHz"),
            ("feedback div", f"{self.fb_div}"),
            ("f_PD / T_PD", f"{self.f_pd/1e9:.4f} GHz / {self.t_pd*1e12:.2f} ps"),
            ("C_SAR", f"{f.c_sar*1e15:.0f} fF"),
            ("SR (diff)", f"{f.slew_rate*1e-9:.3f} mV/ps"),
            ("I_R", f"{f.i_ramp*1e6:.1f} uA"),
            ("DAC", f"{f.dac_bits} b ({f.dac_msb_thermo_bits} b thermo + "
                    f"{f.dac_bits-f.dac_msb_thermo_bits} b bin)"),
            ("DAC LSB", f"{self.dac_lsb_time*1e15:.1f} fs = "
                        f"{self.dac_lsb_volt*1e6:.0f} uV"),
            ("DAC full scale", f"{self.dac_full_scale_volt*1e3:.0f} mV diff "
                               f"({self.dac_full_scale_volt*5e2:.0f} mV s.e.)"),
            ("ADC", f"{f.adc_bits} b, LSB {f.adc_lsb*1e6:.0f} uV "
                    f"= {f.adc_lsb_time*1e15:.1f} fs"),
            ("ADC range", f"{f.adc_full_scale*1e3:.1f} mV = "
                          f"{f.adc_full_scale/f.slew_rate*1e12:.1f} ps"),
            ("effective PD res.", f"{self.effective_bits:.2f} bit / T_PD"),
            ("loop BW", f"{self.loop.bandwidth/1e6:.2f} MHz, zeta="
                        f"{self.loop.damping:.2f}"),
        ]
        w = max(len(k) for k, _ in rows)
        return "\n".join(f"  {k:<{w}} : {v}" for k, v in rows)

    def with_fcw(self, fcw: float) -> "DesignParams":
        return replace(self, fcw=fcw)


def default_design() -> DesignParams:
    """Pure-physics design point: only equations (1)-(13) and published values.

    In-band noise comes out 4 to 7 dB below the measurement.  That is expected
    -- Sec. V says so explicitly: "the main difference between 1 and 100 kHz
    comes from the fact that, in the s-domain analysis, (major) flicker noise
    contribution (from the ramp generator) is neglected."  Use
    :func:`measured_fit_design` to put that term back.
    """
    return DesignParams()


def measured_fit_design() -> DesignParams:
    """Design point that reproduces the measured Fig. 11 / Fig. 12 profiles.

    Fitted by least squares against the four Fig. 11 markers plus both quoted
    jitter figures.  Exactly *one* parameter beyond the published equations is
    needed:

        fdvpd.flicker_corner = 520 kHz

    i.e. the ramp-generator 1/f noise that the paper itself identifies as the
    missing term.  The optimiser was also free to add a white ``excess_noise_db``
    and converged to 0.0 dB, so the published thermal equations need no fudge.

    Resulting agreement (integer-N): 83.5 fs vs 82 fs measured, -115.7 dBc/Hz
    vs -116 dBc/Hz at 110 kHz.  Fractional-N: 97 fs vs 101 fs measured.
    """
    d = DesignParams()
    return replace(d, fdvpd=replace(d.fdvpd, flicker_corner=520e3))


def fractional_design(frac_bit: int = 5, base: DesignParams | None = None
                      ) -> DesignParams:
    """Near-integer fractional channel with only bit ``frac_bit`` asserted.

    The fractional word is applied at the phase detector, so the total FCW
    moves by ``fb_div * 2**-frac_bit``.  ``frac_bit=5`` gives FCW_pd = 20.96875,
    the "FCW = 20.9688 x 2" channel of Fig. 11; sweeping ``frac_bit`` from 2 to
    16 reproduces Fig. 15.
    """
    d = DesignParams() if base is None else base
    return replace(d, fcw=d.fcw - 2.0 ** (-frac_bit) * d.fb_div)
