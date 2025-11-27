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

"""Test suite for utils in pyorbital.orbital."""

import numpy as np
import pytest

from pyorbital.orbital import _get_max_parab, _get_root, compute_azimuth_elevation, ecef_to_topocentric


@pytest.mark.parametrize(
    ("fun", "a", "b", "expected"),
    [
        (lambda x: x**2 - 4, 1, 3, 2.0),  # Basic root test
        (lambda x: x**2 - 4, -3, -1, -2.0),  # Negative root test
        (lambda x: np.sin(x), 3, 4, np.pi),  # High precision near pi
    ],
)
def test_get_root_parametrized(fun, a, b, expected):
    """Test root finding with various functions and intervals."""
    root = _get_root(fun, a, b, tol=1e-6)
    assert np.isclose(root, expected, atol=1e-3)


def test_get_root_unbracketed():
    """Test that root finding fails when no root is bracketed."""

    def fun(x):
        return x**2 + 1

    with pytest.raises(ValueError, match=".*opposite signs.*"):
        _get_root(fun, -1, 1)


@pytest.mark.parametrize(
    ("fun", "a", "b", "expected", "atol"),
    [
        (lambda x: -(x**2) + 4 * x, 1, 3, 2.0, 1e-3),  # Basic peak test
        (lambda x: 5, 0, 10, 5.0, 1.0),  # Flat function
        (lambda x: -np.cos(x), 0, np.pi, np.pi / 2, 1e-2),  # Cosine peak
    ],
)
def test_get_max_parab_parametrized(fun, a, b, expected, atol):
    """Test peak finding with various functions and intervals."""
    peak = _get_max_parab(fun, a, b)
    assert np.isclose(peak, expected, atol=atol)


def test_get_max_parab_divergence():
    """Test fallback behavior for a rising exponential with no peak."""

    def fun(x):
        return np.exp(x)

    peak = _get_max_parab(fun, 0, 1)
    assert 0 <= peak <= 1


def test_ecef_to_topocentric_scalar():
    """Test topocentric conversion with scalar inputs."""
    rx, ry, rz = 1.0, 2.0, 3.0
    lat, lon, theta = np.deg2rad(45), np.deg2rad(60), np.deg2rad(30)
    top_s, top_e, top_z = ecef_to_topocentric(rx, ry, rz, lat, lon, theta)

    assert isinstance(top_s, float)
    assert isinstance(top_e, float)
    assert isinstance(top_z, float)


def test_ecef_to_topocentric_array():
    """Test topocentric conversion with array inputs."""
    rx = np.array([1.0, 2.0])
    ry = np.array([2.0, 3.0])
    rz = np.array([3.0, 4.0])
    lat = np.deg2rad(np.array([45.0, 46.0]))
    lon = np.deg2rad(np.array([60.0, 61.0]))
    theta = np.deg2rad(np.array([30.0, 31.0]))

    top_s, top_e, top_z = ecef_to_topocentric(rx, ry, rz, lat, lon, theta)
    assert top_s.shape == (2,)
    assert top_e.shape == (2,)
    assert top_z.shape == (2,)


def test_compute_azimuth_elevation_scalar():
    """Test azimuth and elevation computation with scalar inputs."""
    top_s, top_e, top_z = 1.0, 1.0, 1.0
    rg = np.sqrt(3)
    az, el = compute_azimuth_elevation(top_s, top_e, top_z, rg)

    assert isinstance(az, float)
    assert isinstance(el, float)
    assert 0 <= az <= 360
    assert -90 <= el <= 90


@pytest.mark.parametrize(
    ("top_s", "top_e", "top_z"),
    [
        (
            np.array([1.0, 0.0]),
            np.array([0.0, 1.0]),
            np.array([1.0, 1.0]),
        ),  # Basic array test
    ],
)
def test_compute_azimuth_elevation_array(top_s, top_e, top_z):
    """Test azimuth and elevation computation with array inputs."""
    rg = np.sqrt(top_s**2 + top_e**2 + top_z**2)
    az, el = compute_azimuth_elevation(top_s, top_e, top_z, rg)
    assert az.shape == (2,)
    assert el.shape == (2,)
    assert np.all((az >= 0) & (az <= 360))
    assert np.all((el >= -90) & (el <= 90))


def test_elevation_clipping():
    """Test elevation clipping when top_z slightly exceeds rg."""
    top_s = np.array([0.0])
    top_e = np.array([0.0])
    top_z = np.array([1.0000001])
    rg = np.array([1.0])

    az, el = compute_azimuth_elevation(top_s, top_e, top_z, rg)
    assert np.isclose(el, 90.0)


def test_zero_distance():
    """Test behavior when observer and satellite are at the same location."""
    rx = ry = rz = 0.0
    rg = 1e-10  # Avoid division by zero
    az, el = compute_azimuth_elevation(rx, ry, rz, rg)
    assert np.isclose(el, 0.0) or np.isclose(el, 90.0)


def test_negative_altitude():
    """Validate azimuth/elevation output when observer is below sea level."""
    rx, ry, rz = 1.0, 1.0, 1.0
    lat, lon, theta = np.deg2rad(0), np.deg2rad(0), np.deg2rad(0)
    top_s, top_e, top_z = ecef_to_topocentric(rx, ry, rz, lat, lon, theta)
    rg = np.sqrt(rx**2 + ry**2 + rz**2)
    az, el = compute_azimuth_elevation(top_s, top_e, top_z, rg)
    assert 0 <= az <= 360
    assert -90 <= el <= 90


def test_polar_observer():
    """Check azimuth/elevation computation near the geographic poles."""
    lat = np.deg2rad(89.999)
    lon = np.deg2rad(0)
    theta = np.deg2rad(0)
    rx, ry, rz = 1.0, 0.0, 0.0
    top_s, top_e, top_z = ecef_to_topocentric(rx, ry, rz, lat, lon, theta)
    rg = np.sqrt(rx**2 + ry**2 + rz**2)
    az, el = compute_azimuth_elevation(top_s, top_e, top_z, rg)
    assert 0 <= az <= 360
    assert -90 <= el <= 90


def test_large_array_input():
    """Stress test azimuth/elevation computation with large input arrays."""
    n = 10000
    rx = np.ones(n)
    ry = np.ones(n)
    rz = np.ones(n)
    lat = np.deg2rad(np.full(n, 45.0))
    lon = np.deg2rad(np.full(n, 60.0))
    theta = np.deg2rad(np.full(n, 30.0))
    top_s, top_e, top_z = ecef_to_topocentric(rx, ry, rz, lat, lon, theta)
    rg = np.sqrt(rx**2 + ry**2 + rz**2)
    az, el = compute_azimuth_elevation(top_s, top_e, top_z, rg)
    assert az.shape == (n,)
    assert el.shape == (n,)


def test_nan_inf_input():
    """Ensure azimuth/elevation handles NaN and Inf inputs gracefully."""
    rx = np.array([np.nan, np.inf])
    ry = np.array([1.0, 1.0])
    rz = np.array([1.0, 1.0])
    rg = np.sqrt(rx**2 + ry**2 + rz**2)
    az, el = compute_azimuth_elevation(rx, ry, rz, rg)
    assert np.isnan(az[0])
    assert np.isfinite(az[1])
    assert np.isnan(el[0])
    assert np.isfinite(el[1])


def test_empty_array_input():
    """Ensure azimuth/elevation handles zero-length arrays gracefully."""
    top_s = np.array([])
    top_e = np.array([])
    top_z = np.array([])
    rg = np.array([])

    az, el = compute_azimuth_elevation(top_s, top_e, top_z, rg)
    assert az.size == 0
    assert el.size == 0


@pytest.mark.parametrize(
    ("lat_deg", "lon_deg"),
    [
        (90.0, 0.0),  # North Pole
        (-90.0, 0.0),  # South Pole
        (0.0, 180.0),  # International Date Line East
        (0.0, -180.0),  # International Date Line West
    ],
)
def test_extreme_lat_lon(lat_deg, lon_deg):
    """Test topocentric conversion at extreme latitude/longitude values."""
    rx, ry, rz = 1.0, 1.0, 1.0
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    theta = 0.0
    top_s, top_e, top_z = ecef_to_topocentric(rx, ry, rz, lat, lon, theta)
    rg = np.sqrt(rx**2 + ry**2 + rz**2)
    az, el = compute_azimuth_elevation(top_s, top_e, top_z, rg)
    assert 0 <= az <= 360
    assert -90 <= el <= 90
