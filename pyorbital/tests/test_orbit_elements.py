#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2012-2024 Pytroll Community

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Test the orbit orbit elements."""

import datetime

import numpy as np
import pytest

from pyorbital.orbital import (
    AE,
    ECC_EPS,
    ECC_LIMIT_HIGH,
    XKMPER,
    XMNPDA,
    OrbitElements,
)


class MockTLE:
    """Mock TLE object for testing OrbitElements."""

    def __init__(self):
        """Initialize mock TLE values for testing."""
        self.epoch = datetime.datetime(2025, 9, 30, 12, 0, 0)
        self.excentricity = 0.001
        self.inclination = 98.7
        self.right_ascension = 120.0
        self.arg_perigee = 87.0
        self.mean_anomaly = 0.0
        self.mean_motion = 14.2
        self.mean_motion_derivative = 0.0
        self.mean_motion_sec_derivative = 0.0
        self.bstar = 0.0001


def test_orbit_elements_computation():
    """Test basic orbital element computations (SMA, Period, Angle Range)."""
    tle = MockTLE()
    orbit = OrbitElements(tle)

    assert isinstance(orbit.semi_major_axis, float)
    assert orbit.semi_major_axis > 0
    assert orbit.period > 0
    assert -np.pi <= orbit.right_ascension_lon <= np.pi


def test_zero_excentricity():
    """Test perigee calculation with zero eccentricity."""
    tle = MockTLE()
    tle.excentricity = 0.0
    orbit = OrbitElements(tle)
    # Perigee == Apogee altitude when e=0
    assert orbit.perigee == pytest.approx(
        (orbit.semi_major_axis / AE - AE) * XKMPER, abs=1e-3
    )


def test_retrograde_orbit_inclination():
    """Test basic inclination check for retrograde orbit."""
    tle = MockTLE()
    tle.inclination = 120.0
    orbit = OrbitElements(tle)
    assert orbit.inclination > np.pi / 2


def test_ra_lon_wrapping():
    """Test longitude wrapping of right ascension with input > 360 deg."""
    tle = MockTLE()
    tle.right_ascension = 400.0  # degrees, > 360
    orbit = OrbitElements(tle)
    assert -np.pi <= orbit.right_ascension_lon <= np.pi


def test_mean_motion_conversion():
    """Test conversion of mean motion to radians per minute."""
    tle = MockTLE()
    orbit = OrbitElements(tle)
    expected = tle.mean_motion * (2 * np.pi / XMNPDA)
    assert orbit.mean_motion == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize("incl_deg", [0, 45, 90, 135])
def test_semi_major_axis_vs_inclination(incl_deg):
    """Test semi-major axis computation across inclinations."""
    tle = MockTLE()
    tle.inclination = incl_deg
    orbit = OrbitElements(tle)
    assert orbit.semi_major_axis > 0


@pytest.mark.parametrize("bstar", [0.0, 1e-5, 0.0001, 0.01])
def test_bstar_scaling(bstar):
    """Test scaling of bstar drag term."""
    tle = MockTLE()
    tle.bstar = bstar
    orbit = OrbitElements(tle)
    assert orbit.bstar == pytest.approx(bstar * AE)


def test_mean_motion_derivatives():
    """Test mean motion derivatives are correctly computed."""
    tle = MockTLE()
    orbit = OrbitElements(tle)
    tle.mean_motion_derivative = 1e-5
    tle.mean_motion_sec_derivative = 1e-7
    orbit = OrbitElements(tle)
    assert orbit.mean_motion_derivative > 0
    assert orbit.mean_motion_sec_derivative > 0


def test_apogee_computation():
    """Test apogee altitude calculation."""
    tle = MockTLE()
    orbit = OrbitElements(tle)
    expected_apogee = (
        (orbit.semi_major_axis * (1 + orbit.excentricity)) / AE - AE
    ) * XKMPER
    assert orbit.apogee == pytest.approx(expected_apogee, abs=1e-3)


def test_semi_major_axis_accuracy():
    """Test semi-major axis against analytical reference value (in km)."""
    tle = MockTLE()
    orbit = OrbitElements(tle)
    expected_sma_km = 7200.645659667062
    computed_sma_km = orbit.semi_major_axis * XKMPER / AE
    assert computed_sma_km == pytest.approx(expected_sma_km, abs=0.1)


def test_orbital_period_accuracy():
    """Test orbital period against corrected mean motion (in seconds)."""
    tle = MockTLE()
    orbit = OrbitElements(tle)
    expected_period_sec = 6080.901176943267
    computed_period_sec = orbit.period * 60
    assert computed_period_sec == pytest.approx(expected_period_sec, abs=1.0)


def test_velocity_at_perigee_accuracy():
    """Test orbital velocity at perigee against analytical reference."""
    tle = MockTLE()
    orbit = OrbitElements(tle)
    expected_velocity_kms = 7.44762253625217
    computed_velocity_kms = orbit.velocity_at_perigee()
    assert computed_velocity_kms == pytest.approx(expected_velocity_kms, abs=0.005)


def test_velocity_at_apogee_accuracy():
    """Test orbital velocity at apogee against analytical reference."""
    tle = MockTLE()
    orbit = OrbitElements(tle)
    expected_velocity_kms = 7.432742171544375
    computed_velocity_kms = orbit.velocity_at_apogee()
    assert computed_velocity_kms == pytest.approx(expected_velocity_kms, abs=0.005)


def test_position_vector_in_orbital_plane_perigee_accuracy():
    """Test position vector at perigee (Mean Anomaly = 0°)."""
    tle = MockTLE()
    tle.excentricity = 0.001
    tle.mean_anomaly = 0.0
    orbit = OrbitElements(tle)
    pos = orbit.position_vector_in_orbital_plane()
    expected_x = 1.1278289051591721  # Earth radii
    expected_y = 0.0
    assert pos[0] == pytest.approx(expected_x, rel=1e-6)
    assert pos[1] == pytest.approx(expected_y, abs=1e-8)


@pytest.mark.parametrize(
    ("mean_anomaly_deg", "expected_radius"),
    [
        (0, 7193.445014007397),
        (90, 7200.6528603079205),
        (180, 7207.84630532673),
        (270, 7200.6528603079205),
    ],
)
def test_position_vector_in_orbital_plane_varied(mean_anomaly_deg, expected_radius):
    """Verify position vector magnitude matches expected radius at given mean anomaly."""
    tle = MockTLE()
    tle.mean_anomaly = mean_anomaly_deg
    orbit = OrbitElements(tle)
    pos = orbit.position_vector_in_orbital_plane()
    actual_radius = np.linalg.norm(pos) * XKMPER
    assert np.isclose(actual_radius, expected_radius, rtol=1e-6)


@pytest.mark.parametrize(
    ("excentricity", "expected"),
    [
        (0.0005, True),
        (0.01, False),
    ],
)
def test_is_circular_property(excentricity, expected):
    """Test circular orbit detection based on eccentricity."""
    tle = MockTLE()
    tle.excentricity = excentricity
    orbit = OrbitElements(tle)
    assert orbit.is_circular == expected


@pytest.mark.parametrize(
    ("inclination_deg", "expected"),
    [
        (100.0, True),
        (80.0, False),
    ],
)
def test_is_retrograde_property(inclination_deg, expected):
    """Test retrograde orbit detection based on inclination."""
    tle = MockTLE()
    tle.inclination = inclination_deg
    orbit = OrbitElements(tle)
    assert orbit.is_retrograde == expected


@pytest.mark.parametrize("excentricity", [0.001, 0.01, 0.1, 0.5])
def test_velocity_perigee_greater_than_apogee(excentricity):
    """Ensure velocity at perigee is greater than at apogee for elliptical orbits."""
    tle = MockTLE()
    tle.excentricity = excentricity
    orbit = OrbitElements(tle)
    v_perigee = orbit.velocity_at_perigee()
    v_apogee = orbit.velocity_at_apogee()
    assert v_perigee > v_apogee


@pytest.mark.parametrize("excentricity", [0.0, ECC_EPS / 10, ECC_EPS])
def test_velocity_equal_for_circular_orbits(excentricity):
    """Ensure velocity at perigee equals velocity at apogee for nearly circular orbits."""
    tle = MockTLE()
    tle.excentricity = excentricity
    orbit = OrbitElements(tle)
    v_perigee = orbit.velocity_at_perigee()
    v_apogee = orbit.velocity_at_apogee()
    assert v_perigee == pytest.approx(v_apogee, rel=1e-5)


def test_high_eccentricity_limit():
    """Test behavior near the upper eccentricity limit."""
    tle = MockTLE()
    tle.excentricity = ECC_LIMIT_HIGH
    orbit = OrbitElements(tle)
    assert orbit.semi_major_axis > 0
    assert orbit.velocity_at_perigee() > orbit.velocity_at_apogee()


def test_sgp4_unnormalization_value():
    """Validate SGP4 un-normalization of mean motion and semi-major axis.

    This test uses hardcoded reference values from validated SGP4 output.
    """

    class RefTLE(MockTLE):
        """Reference TLE with fixed orbital parameters for precision testing."""

        def __init__(self):
            super().__init__()
            """Initialize reference TLE values for SGP4 unnormalization test."""
            self.inclination = 98.7408  # degrees
            self.excentricity = 0.001140
            self.mean_motion = 14.28634842  # rev/day
            self.epoch = datetime.datetime(2200, 1, 1, 0, 0, 0)

    tle = RefTLE()
    orbit = OrbitElements(tle)

    # Hardcoded expected values from validated SGP4 propagation
    expected_mean_motion = 0.06237319306246916  # rad/min
    expected_semi_major_axis = 1.1244009310620886  # Earth radii

    assert orbit.original_mean_motion == pytest.approx(expected_mean_motion, rel=1e-8)
    assert orbit.semi_major_axis == pytest.approx(expected_semi_major_axis, rel=1e-8)
