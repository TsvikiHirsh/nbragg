"""Shared pytest configuration for the nbragg test suite.

Several tests build small synthetic datasets with NumPy's legacy global RNG
(``np.random.randint`` / ``np.random.random`` etc.) without seeding it. Left
unseeded, the random data occasionally lands on a degenerate case where a fit
does not converge, making those tests flaky across NumPy/SciPy/lmfit/NCrystal
versions (e.g. between the pinned local environment and the latest packages CI
installs). Seeding the global RNG before every test makes the whole suite
deterministic without touching each individual test.
"""
import numpy as np
import pytest


@pytest.fixture(autouse=True)
def _seed_global_rng():
    """Seed NumPy's legacy global RNG before each test for reproducibility."""
    np.random.seed(0)
    yield
