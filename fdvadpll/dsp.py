"""Noise synthesis and spectral-analysis helpers.

Conventions used throughout the package
---------------------------------------
* ``S_x(f)``  one-sided power spectral density, ``f`` in [0, fs/2].
  The variance of a real sequence is ``integral S_x df``.
* ``L(f)``    single-sideband phase noise, ``L = S_phi / 2`` with ``S_phi``
  the one-sided PSD of the excess phase in radians.
* A sequence of uncorrelated period deviations with standard deviation
  ``sigma_dT`` at carrier ``f0`` produces ``L(f) = f0**3 * sigma_dT**2 / f**2``.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "colored_noise",
    "one_over_f",
    "time_error_from_phase_noise",
    "dco_period_jitter",
    "welch_psd",
    "phase_noise_from_time_error",
    "integrate_phase_noise",
    "rms_jitter",
    "find_spurs",
    "spur_level",
]


# --------------------------------------------------------------------------
#  Synthesis
# --------------------------------------------------------------------------
def colored_noise(n: int, fs: float, psd_func, rng: np.random.Generator,
                  f_min: float | None = None) -> np.ndarray:
    """Real sequence of length ``n`` whose one-sided PSD is ``psd_func(f)``.

    Built by shaping white noise in the frequency domain.  ``psd_func`` is
    evaluated on the rFFT grid; the DC bin is forced to zero and everything
    below ``f_min`` (default: the first non-zero bin) is flattened so that
    divergent laws such as 1/f**3 stay finite.
    """
    n = int(n)
    nf = n // 2 + 1
    f = np.arange(nf) * fs / n
    df = fs / n
    f_min = df if f_min is None else max(f_min, df)
    f_eval = np.maximum(f, f_min)
    s = np.asarray(psd_func(f_eval), dtype=float)
    s[0] = 0.0

    # One-sided PSD of an rFFT: S(f_k) = 2 |X_k|^2 / (n*fs)
    mag = np.sqrt(s * n * fs / 2.0)
    phase = rng.uniform(0.0, 2.0 * np.pi, nf)
    x = mag * np.exp(1j * phase)
    x[0] = 0.0
    if n % 2 == 0:
        # Nyquist bin must be real; it also carries half the two-sided power.
        x[-1] = mag[-1] / math.sqrt(2.0) * np.sign(rng.standard_normal())
    return np.fft.irfft(x, n)


def one_over_f(n: int, fs: float, amplitude: float,
               rng: np.random.Generator, f_min: float | None = None):
    """Sequence with one-sided PSD ``amplitude / f``."""
    return colored_noise(n, fs, lambda f: amplitude / f, rng, f_min)


def time_error_from_phase_noise(n: int, fs: float, f0: float, pn_func,
                                rng: np.random.Generator) -> np.ndarray:
    """Timing-error sequence [s] realising a given ``L(f)`` (linear) profile.

    ``S_x(f) = 2 L(f) / (2*pi*f0)**2``
    """
    scale = 2.0 / (2.0 * math.pi * f0) ** 2
    return colored_noise(n, fs, lambda f: scale * np.asarray(pn_func(f)), rng)


def dco_period_jitter(n_edges: int, f_edge: float, sigma_white: float,
                      f_corner: float, rng: np.random.Generator) -> np.ndarray:
    """Per-edge period deviations [s] giving 1/f**2 + 1/f**3 phase noise.

    White period jitter of ``sigma_white`` yields ``L = f0**3 sigma**2 / f**2``.
    For the flicker part to cross it at ``f_corner`` the period-deviation PSD
    must be ``S_dT(f) = 2 sigma**2 f_corner / (f0 * f)``  (see module docstring
    for the L <-> S_dT mapping).
    """
    white = rng.standard_normal(n_edges) * sigma_white
    if f_corner <= 0.0:
        return white
    amp = 2.0 * sigma_white ** 2 * f_corner / f_edge
    return white + one_over_f(n_edges, f_edge, amp, rng)


# --------------------------------------------------------------------------
#  Measurement
# --------------------------------------------------------------------------
def welch_psd(x: np.ndarray, fs: float, nperseg: int | None = None,
              noverlap: int | None = None, window: str = "hann",
              detrend: bool = True):
    """One-sided PSD via Welch's method.  Returns ``(f, S)`` excluding DC."""
    from scipy import signal

    x = np.asarray(x, dtype=float)
    if nperseg is None:
        nperseg = min(len(x), 1 << int(math.log2(max(len(x) // 8, 256))))
    f, s = signal.welch(x, fs=fs, nperseg=nperseg, noverlap=noverlap,
                        window=window, detrend="constant" if detrend else False,
                        return_onesided=True, scaling="density")
    return f[1:], s[1:]


def phase_noise_from_time_error(x: np.ndarray, fs: float, f0: float,
                                nperseg: int | None = None):
    """``L(f)`` [linear] from a timing-error sequence sampled at ``fs``.

    ``S_phi = (2*pi*f0)**2 * S_x``,  ``L = S_phi/2``.
    """
    f, s_x = welch_psd(x, fs, nperseg=nperseg)
    s_phi = (2.0 * math.pi * f0) ** 2 * s_x
    return f, s_phi / 2.0


def integrate_phase_noise(f, L_lin, f_lo: float, f_hi: float) -> float:
    """Integrated phase noise ``2*int L df`` [rad^2] over ``f_lo .. f_hi``."""
    f = np.asarray(f, dtype=float)
    L_lin = np.asarray(L_lin, dtype=float)
    m = (f >= f_lo) & (f <= f_hi)
    if m.sum() < 2:
        raise ValueError(f"band {f_lo:g}..{f_hi:g} Hz has < 2 spectral points "
                         f"(resolution {f[1]-f[0]:g} Hz)")
    return 2.0 * float(np.trapezoid(L_lin[m], f[m]))


def rms_jitter(f, L_lin, f0: float, f_lo: float = 10e3,
               f_hi: float = 40e6) -> float:
    """RMS jitter [s] over the given integration band."""
    ipn = integrate_phase_noise(f, L_lin, f_lo, f_hi)
    return math.sqrt(ipn) / (2.0 * math.pi * f0)


def find_spurs(x: np.ndarray, fs: float, f0: float, n_fft: int | None = None,
               exclude_dc_bins: int = 4):
    """Periodogram of a timing-error sequence, in dBc.

    Returns ``(f, dbc)`` where ``dbc`` is the power in each bin relative to the
    carrier, i.e. what a spectrum analyser reads for a discrete spur.  Uses a
    Blackman-Harris window and corrects for its coherent gain so that a pure
    tone reads at its true level.
    """
    from scipy.signal.windows import blackmanharris

    x = np.asarray(x, dtype=float)
    n = len(x) if n_fft is None else int(n_fft)
    x = x[:n] - x[:n].mean()
    w = blackmanharris(n)
    cg = w.sum() / n                       # coherent gain
    phi = 2.0 * math.pi * f0 * x           # excess phase, rad
    spec = np.fft.rfft(phi * w) / (n * cg)
    # A real tone of peak phase deviation A lands in one rFFT bin with
    # |spec| = A/2, which is already the sideband-to-carrier amplitude ratio of
    # narrow-band PM (J_1(A)/J_0(A) -> A/2), i.e. what an analyser reads.
    p = np.abs(spec) ** 2
    f = np.arange(len(spec)) * fs / n
    p[:exclude_dc_bins] = 0.0
    with np.errstate(divide="ignore"):
        return f[1:], 10.0 * np.log10(p[1:])


def spur_level(x: np.ndarray, fs: float, f0: float, f_spur: float,
               search_bins: int = 3, n_fft: int | None = None):
    """Level [dBc] of the spur nearest ``f_spur``.  Returns ``(f_found, dbc)``."""
    f, dbc = find_spurs(x, fs, f0, n_fft)
    idx = int(np.argmin(np.abs(f - f_spur)))
    lo = max(idx - search_bins, 0)
    hi = min(idx + search_bins + 1, len(f))
    j = lo + int(np.argmax(dbc[lo:hi]))
    return float(f[j]), float(dbc[j])
