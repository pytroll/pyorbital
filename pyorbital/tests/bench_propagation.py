"""Benchmarks for SGP4/SDP4 orbit propagation.

Run with::

    pytest pyorbital/tests/bench_propagation.py --benchmark-only -v

Record a baseline before changing the propagator::

    pytest pyorbital/tests/bench_propagation.py --benchmark-only --benchmark-save=pre-sdp4

and compare later work against it::

    pytest pyorbital/tests/bench_propagation.py --benchmark-only \
        --benchmark-compare=pre-sdp4 --benchmark-compare-fail=median:8%

Three near-earth scenarios establish the reference cost of the existing SGP4
path: a single propagation (per-call overhead), a day of propagations as one
array (vectorised throughput), and construction of an ``Orbital`` object (TLE
parsing plus the one-off ``_SGDP4Base`` coefficient computation).
"""

import datetime as dt

import numpy as np
import pytest

# FY-3F, epoch 2026-04-01. Near-earth (~100 min period), the common LEO case.
_TLE1 = "1 57490U 23111A   26091.51900704  .00000079  00000+0  57513-4 0  9997"
_TLE2 = "2 57490  98.6883 162.9408 0000839  82.8428 277.2844 14.19924724137986"
_T0 = dt.datetime(2026, 4, 1, 12, 0, 0)

_VECTOR_POINTS = 10000


@pytest.fixture(scope="module")
def near_earth_orbital():
    """Pre-built Orbital object, so construction cost stays out of the timing."""
    from pyorbital.orbital import Orbital
    return Orbital("FY-3F", line1=_TLE1, line2=_TLE2)


@pytest.fixture(scope="module")
def one_day_of_times():
    """A day of propagation times, evenly spaced, as a datetime64 array."""
    start = np.datetime64(_T0)
    offsets = np.linspace(0, 86400, _VECTOR_POINTS)
    return start + (offsets * 1e6).astype("timedelta64[us]")


def _construct_orbital():
    from pyorbital.orbital import Orbital
    return Orbital("FY-3F", line1=_TLE1, line2=_TLE2)


def _assert_plausible_positions(pos):
    """Guard against silently timing a wrong result."""
    radius = np.sqrt(np.sum(np.asarray(pos) ** 2, axis=0))
    assert np.all(np.isfinite(radius))
    assert np.all(radius > 6378.0)
    assert np.all(radius < 7500.0)


def test_bench_near_earth_scalar(benchmark, near_earth_orbital):
    """Benchmark: single near-earth propagation (per-call overhead)."""
    pos, vel = benchmark(near_earth_orbital.get_position, _T0, False)

    _assert_plausible_positions(pos)


def test_bench_near_earth_vector_day(benchmark, near_earth_orbital, one_day_of_times):
    """Benchmark: a day of near-earth propagations as one array."""
    pos, vel = benchmark(near_earth_orbital.get_position, one_day_of_times, False)

    assert np.asarray(pos).shape == (3, _VECTOR_POINTS)
    _assert_plausible_positions(pos)


def test_bench_orbital_init(benchmark):
    """Benchmark: Orbital construction (TLE parsing and orbit coefficients)."""
    orbital = benchmark(_construct_orbital)

    pos, vel = orbital.get_position(_T0, False)
    _assert_plausible_positions(pos)
