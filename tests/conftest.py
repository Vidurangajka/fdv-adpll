"""Shared fixtures.

Simulation-based tests are expensive, so the few runs that more than one test
wants to look at are built once per session and cached.
"""

import numpy as np
import pytest

from fdvadpll import DesignParams, FdvPll, default_design, fractional_design


@pytest.fixture(scope="session")
def design() -> DesignParams:
    return default_design()


@pytest.fixture(scope="session")
def rng() -> np.random.Generator:
    return np.random.default_rng(12345)


@pytest.fixture(scope="session")
def integer_run():
    """Integer-N lock, long enough for the 10 kHz integration limit."""
    return FdvPll(default_design(), seed=7).run(1 << 16)


@pytest.fixture(scope="session")
def fractional_run():
    """FCW_pd = 20.96875 -- the near-integer channel of Fig. 11."""
    return FdvPll(fractional_design(5), seed=7).run(1 << 16)
