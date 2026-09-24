"""Behavioural model of a fractional-N DPLL whose phase error is quantised in
the fully differential voltage domain, plus an all-digital (ADPLL) extension.

Reference implementation of

    L. Wu, T. Burger, P. Schoenle and Q. Huang, "A Power-Efficient Fractional-N
    DPLL With Phase Error Quantized in Fully Differential-Voltage Domain,"
    IEEE J. Solid-State Circuits, vol. 56, no. 4, pp. 1254-1264, April 2021.

Quick start
-----------
>>> from fdvadpll import DesignParams, FdvPll
>>> pll = FdvPll(DesignParams(), seed=1)
>>> res = pll.run(1 << 15)
>>> print(res.summary())            # doctest: +SKIP
"""

from .params import (DesignParams, DcoParams, FdvpdParams, LoopParams,
                     PowerParams, RefParams, default_design, fractional_design,
                     measured_fit_design)
from .blocks import CurrentDAC, Fdvpd, RampGenerator, SarAdc
from .dco import Dco
from .dlf import DigitalLoopFilter, FrequencyLockLoop
from .pll import FdvPll, SimResult
from .sdm import MashSdm
from . import calib, dsp, noise

__version__ = "1.0.0"

__all__ = [
    "DesignParams", "DcoParams", "FdvpdParams", "LoopParams", "PowerParams",
    "RefParams", "default_design", "fractional_design", "measured_fit_design",
    "CurrentDAC", "Fdvpd", "RampGenerator", "SarAdc",
    "Dco", "DigitalLoopFilter", "FrequencyLockLoop",
    "FdvPll", "SimResult", "MashSdm",
    "calib", "dsp", "noise",
]
