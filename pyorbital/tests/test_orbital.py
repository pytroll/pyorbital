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

"""Test the geoloc orbital."""

import datetime as dt
import unittest
import warnings
from unittest import mock

import dask.array as da
import numpy as np
import pytest
import xarray as xr

from pyorbital.orbital import Orbital, get_observer_look

eps_deg = 10e-3


class Test(unittest.TestCase):
    """Basic test class for unittesting the pyorbital.orbital class."""

    def test_get_orbit_number(self):
        """Testing getting the orbitnumber from the TLEs."""
        sat = Orbital(
            "NPP",
            line1="1 37849U 11061A   12017.90990040 " "-.00000112  00000-0 -32693-4 0   772",
            line2="2 37849  98.7026 317.8811 0001845  " "92.4533 267.6830 14.19582686 11574",
        )
        dobj = dt.datetime(2012, 1, 18, 8, 4, 19)
        orbnum = sat.get_orbit_number(dobj)
        assert orbnum == 1163

    def test_sublonlat(self):
        """Test getting the sub-satellite position."""
        sat = Orbital(
            "ISS (ZARYA)",
            line1="1 25544U 98067A   03097.78853147  " ".00021906  00000-0  28403-3 0  8652",
            line2="2 25544  51.6361  13.7980 0004256  " "35.6671  59.2566 15.58778559250029",
        )
        d = dt.datetime(2003, 3, 23, 0, 3, 22)
        lon, lat, alt = sat.get_lonlatalt(d)
        expected_lon = -68.199894472013213
        expected_lat = 23.159747677881075
        expected_alt = 392.01953430856935
        assert np.abs(lon - expected_lon) < eps_deg, "Calculation of sublon failed"
        assert np.abs(lat - expected_lat) < eps_deg, "Calculation of sublat failed"
        assert np.abs(alt - expected_alt) < eps_deg, "Calculation of altitude failed"

    def test_observer_look(self):
        """Test getting the observer look angles."""
        sat = Orbital(
            "ISS (ZARYA)",
            line1="1 25544U 98067A   03097.78853147  " ".00021906  00000-0  28403-3 0  8652",
            line2="2 25544  51.6361  13.7980 0004256  " "35.6671  59.2566 15.58778559250029",
        )
        d = dt.datetime(2003, 3, 23, 0, 3, 22)
        az, el = sat.get_observer_look(d, -84.39733, 33.775867, 0)
        expected_az = 122.45169655331965
        expected_el = 1.9800219611255456
        assert np.abs(az - expected_az) < eps_deg, "Calculation of azimut failed"
        assert np.abs(el - expected_el) < eps_deg, "Calculation of elevation failed"

    def test_orbit_num_an(self):
        """Test getting orbit number - ascending node."""
        sat = Orbital(
            "METOP-A",
            line1="1 29499U 06044A   11254.96536486  " ".00000092  00000-0  62081-4 0  5221",
            line2="2 29499  98.6804 312.6735 0001758 " "111.9178 248.2152 14.21501774254058",
        )
        d = dt.datetime(2011, 9, 14, 5, 30)
        assert sat.get_orbit_number(d) == 25437

    def test_orbit_num_non_an(self):
        """Test getting orbit number - not ascending node."""
        sat = Orbital(
            "METOP-A",
            line1="1 29499U 06044A   13060.48822809  " ".00000017  00000-0  27793-4 0  9819",
            line2="2 29499  98.6639 121.6164 0001449  " "71.9056  43.3132 14.21510544330271",
        )
        dt = np.timedelta64(98, "m")
        assert sat.get_orbit_number(sat.tle.epoch + dt) == 33028

    def test_orbit_num_equator(self):
        """Test getting orbit numbers when being around equator."""
        sat = Orbital(
            "SUOMI NPP",
            line1="1 37849U 11061A   13061.24611272  " ".00000048  00000-0  43679-4 0  4334",
            line2="2 37849  98.7444   1.0588 0001264  " "63.8791 102.8546 14.19528338 69643",
        )
        t1 = dt.datetime(2013, 3, 2, 22, 2, 25)
        t2 = dt.datetime(2013, 3, 2, 22, 3, 00)
        on1 = sat.get_orbit_number(t1)
        on2 = sat.get_orbit_number(t2)
        assert on2 >= on1
        assert on2 - on1 <= 1
        pos1, vel1 = sat.get_position(t1, normalize=False)
        pos2, vel2 = sat.get_position(t2, normalize=False)
        del vel1, vel2
        assert pos1[2] < 0
        assert pos2[2] > 0

    def test_get_next_passes_apogee(self):
        """Regression test #22."""
        line1 = "1 24793U 97020B   18065.48735489  " ".00000075  00000-0  19863-4 0  9994"
        line2 = "2 24793  86.3994 209.3241 0002020  " "89.8714 270.2713 14.34246429 90794"

        orb = Orbital("IRIDIUM 7 [+]", line1=line1, line2=line2)
        d = dt.datetime(2018, 3, 7, 3, 30, 15)
        res = orb.get_next_passes(d, 1, 170.556, -43.368, 0.5, horizon=40)
        assert abs(res[0][2] - dt.datetime(2018, 3, 7, 3, 48, 13, 178439)) < dt.timedelta(seconds=0.01)

    def test_get_next_passes_tricky(self):
        """Check issue #34 for reference."""
        line1 = "1 43125U 18004Q   18251.42128650 " "+.00001666 +00000-0 +73564-4 0  9991"

        line2 = "2 43125 097.5269 314.3317 0010735 " "157.6344 202.5362 15.23132245036381"

        orb = Orbital("LEMUR-2-BROWNCOW", line1=line1, line2=line2)
        d = dt.datetime(2018, 9, 8)

        res = orb.get_next_passes(d, 72, -8.174163, 51.953319, 0.05, horizon=5)

        assert abs(res[0][2] - dt.datetime(2018, 9, 8, 9, 5, 46, 375248)) < dt.timedelta(seconds=0.01)
        assert abs(res[-1][2] - dt.datetime(2018, 9, 10, 22, 15, 3, 143469)) < dt.timedelta(seconds=0.01)

        assert len(res) == 15

    def test_get_next_passes_issue_22(self):
        """Check that max."""
        line1 = "1 28654U 05018A   21083.16603416  .00000102  00000-0  79268-4 0  9999"
        line2 = "2 28654  99.0035 147.6583 0014816 159.4931 200.6838 14.12591533816498"

        orb = Orbital("NOAA 18", line1=line1, line2=line2)
        t = dt.datetime(2021, 3, 9, 22)
        next_passes = orb.get_next_passes(t, 1, -15.6335, 27.762, 0.0)
        rise, fall, max_elevation = next_passes[0]
        assert rise < max_elevation < fall

    @mock.patch("pyorbital.orbital.Orbital.get_lonlatalt")
    def test_utc2local(self, get_lonlatalt):
        """Test converting UTC to local time."""
        get_lonlatalt.return_value = -45, None, None
        sat = Orbital(
            "METOP-A",
            line1="1 29499U 06044A   13060.48822809  " ".00000017  00000-0  27793-4 0  9819",
            line2="2 29499  98.6639 121.6164 0001449  " "71.9056  43.3132 14.21510544330271",
        )
        assert sat.utc2local(dt.datetime(2009, 7, 1, 12)) == dt.datetime(2009, 7, 1, 9)

    @mock.patch("pyorbital.orbital.Orbital.utc2local")
    @mock.patch("pyorbital.orbital.Orbital.get_orbit_number")
    def test_get_equatorial_crossing_time(self, get_orbit_number, utc2local):
        """Test get the equatorial crossing time."""

        def get_orbit_number_patched(utc_time, **kwargs):
            utc_time = np.datetime64(utc_time)
            diff = (utc_time - np.datetime64("2009-07-01 12:38:12")) / np.timedelta64(7200, "s")
            return 1234 + diff

        get_orbit_number.side_effect = get_orbit_number_patched
        utc2local.return_value = "local_time"
        sat = Orbital(
            "METOP-A",
            line1="1 29499U 06044A   13060.48822809  " ".00000017  00000-0  27793-4 0  9819",
            line2="2 29499  98.6639 121.6164 0001449  " "71.9056  43.3132 14.21510544330271",
        )

        # Ascending node
        res = sat.get_equatorial_crossing_time(tstart=dt.datetime(2009, 7, 1, 12), tend=dt.datetime(2009, 7, 1, 13))
        exp = dt.datetime(2009, 7, 1, 12, 38, 12)
        assert res - exp < dt.timedelta(seconds=0.01)

        # Descending node
        res = sat.get_equatorial_crossing_time(
            tstart=dt.datetime(2009, 7, 1, 12),
            tend=dt.datetime(2009, 7, 1, 14, 0),
            node="descending",
        )
        exp = dt.datetime(2009, 7, 1, 13, 38, 12)
        assert res - exp < dt.timedelta(seconds=0.01)

        # Conversion to local time
        res = sat.get_equatorial_crossing_time(
            tstart=dt.datetime(2009, 7, 1, 12),
            tend=dt.datetime(2009, 7, 1, 14),
            local_time=True,
        )
        assert res == "local_time"

    def test_get_next_passes_basic(self):
        """Test basic pass detection with known TLE and observer."""
        line1 = "1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652"
        line2 = "2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029"
        sat = Orbital("ISS (ZARYA)", line1=line1, line2=line2)

        observer_lon = -84.39733
        observer_lat = 33.775867
        observer_alt = 0.0
        start_time = dt.datetime(2003, 3, 23, 0, 0)

        passes = sat.get_next_passes(start_time, length=2, lon=observer_lon, lat=observer_lat, alt=observer_alt)
        assert len(passes) > 0
        for rise, fall, peak in passes:
            assert rise < peak < fall

    def test_get_next_passes_no_passes(self):
        """Test when no passes occur within the time window."""
        line1 = "1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652"
        line2 = "2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029"
        sat = Orbital("ISS (ZARYA)", line1=line1, line2=line2)

        observer_lon = 0.0
        observer_lat = -90.0  # South Pole, unlikely for ISS
        observer_alt = 0.0
        start_time = dt.datetime(2003, 3, 23, 0, 0)

        passes = sat.get_next_passes(start_time, length=1, lon=observer_lon, lat=observer_lat, alt=observer_alt)
        assert passes == []

    def test_get_next_passes_with_horizon(self):
        """Test pass detection with a non-zero horizon elevation."""
        line1 = "1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652"
        line2 = "2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029"
        sat = Orbital("ISS (ZARYA)", line1=line1, line2=line2)

        observer_lon = -84.39733
        observer_lat = 33.775867
        observer_alt = 0.0
        start_time = dt.datetime(2003, 3, 23, 0, 0)

        passes = sat.get_next_passes(
            start_time,
            length=2,
            lon=observer_lon,
            lat=observer_lat,
            alt=observer_alt,
            horizon=20,
        )
        for rise, fall, peak in passes:
            assert rise < peak < fall

    def test_get_next_passes_tolerance_effect(self):
        """Test that changing tolerance affects timing precision."""
        line1 = "1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652"
        line2 = "2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029"
        sat = Orbital("ISS (ZARYA)", line1=line1, line2=line2)

        observer_lon = -84.39733
        observer_lat = 33.775867
        observer_alt = 0.0
        start_time = dt.datetime(2003, 3, 23, 0, 0)

        passes_lo_tol = sat.get_next_passes(
            start_time,
            length=2,
            lon=observer_lon,
            lat=observer_lat,
            alt=observer_alt,
            tol=0.001,
        )
        passes_hi_tol = sat.get_next_passes(
            start_time,
            length=2,
            lon=observer_lon,
            lat=observer_lat,
            alt=observer_alt,
            tol=1.0,
        )

        assert len(passes_lo_tol) == len(passes_hi_tol)
        for (r1, f1, p1), (r2, f2, p2) in zip(passes_lo_tol, passes_hi_tol):
            assert abs((r1 - r2).total_seconds()) < 60  # Should be close, but not identical

    def test_observer_look_radian_bug(self):
        """Detect bug caused by passing observer coordinates in radians instead of degrees."""
        line1 = "1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652"
        line2 = "2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029"
        sat = Orbital("ISS (ZARYA)", line1=line1, line2=line2)

        time = dt.datetime(2003, 3, 23, 0, 3, 22)
        obs_lon_deg, obs_lat_deg = -84.39733, 33.775867
        obs_lon_rad, obs_lat_rad = np.deg2rad(obs_lon_deg), np.deg2rad(obs_lat_deg)

        # Correct usage: degrees
        az_deg, el_deg = sat.get_observer_look(time, obs_lon_deg, obs_lat_deg, 0)

        # Incorrect usage: radians
        az_rad, el_rad = sat.get_observer_look(time, obs_lon_rad, obs_lat_rad, 0)

        assert abs(az_deg - az_rad) > 10, "Azimuth unexpectedly similar—possible bug hidden"
        assert abs(el_deg - el_rad) > 10, "Elevation unexpectedly similar—possible bug hidden"

    def test_get_next_passes_radian_bug(self):
        """Detect bug caused by passing observer coordinates in radians instead of degrees."""
        line1 = "1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652"
        line2 = "2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029"
        sat = Orbital("ISS (ZARYA)", line1=line1, line2=line2)

        start_time = dt.datetime(2003, 3, 23, 0, 0)
        observer_lon_deg, observer_lat_deg = -84.39733, 33.775867
        observer_lon_rad, observer_lat_rad = np.deg2rad(observer_lon_deg), np.deg2rad(observer_lat_deg)

        # Correct usage: degrees
        passes_deg = sat.get_next_passes(start_time, length=2, lon=observer_lon_deg, lat=observer_lat_deg, alt=0)

        # Incorrect usage: radians
        passes_rad = sat.get_next_passes(start_time, length=2, lon=observer_lon_rad, lat=observer_lat_rad, alt=0)

        assert len(passes_deg) > 0, "Expected passes not found with correct input"
        assert len(passes_rad) == 0 or any(
            abs((p1[0] - p2[0]).total_seconds()) > 600 for p1, p2 in zip(passes_deg, passes_rad)
        ), "Radian input unexpectedly produced similar results—possible bug hidden"

    def test_find_aos_basic(self):
        """Test basic AOS detection for a known satellite and location."""
        sat = Orbital(
            "ISS (ZARYA)",
            line1="1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
            line2="2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029",
        )
        utc_time = dt.datetime(2003, 3, 23, 0, 0, 0)
        aos = sat.find_aos(utc_time, lon=-84.39733, lat=33.775867, alt=0)
        assert aos is not None, "AOS should be detected"
        assert aos > utc_time, "AOS must be after the input time"

    def test_find_aol_basic(self):
        """Test basic AOL detection for a known satellite and location."""
        sat = Orbital(
            "ISS (ZARYA)",
            line1="1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
            line2="2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029",
        )
        utc_time = dt.datetime(2003, 3, 23, 0, 0, 0)
        aol = sat.find_aol(utc_time, lon=-84.39733, lat=33.775867, alt=0)
        assert aol is not None, "AOL should be detected"
        assert aol > utc_time, "AOL must be after the input time"

    def test_find_aos_no_pass(self):
        """Test AOS detection when no pass occurs within search window."""
        sat = Orbital(
            "ISS (ZARYA)",
            line1="1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
            line2="2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029",
        )
        # Use a location far from satellite's inclination
        utc_time = dt.datetime(2003, 3, 23, 0, 0, 0)
        aos = sat.find_aos(utc_time, lon=0, lat=89.0, alt=0, horizon=10)
        assert aos is None, "No AOS should be found near the poles for ISS"

    def test_find_aos_edge_horizon(self):
        """Test AOS detection with a high horizon angle."""
        sat = Orbital(
            "ISS (ZARYA)",
            line1="1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
            line2="2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029",
        )
        utc_time = dt.datetime(2003, 3, 23, 0, 0, 0)
        aos = sat.find_aos(utc_time, lon=-84.39733, lat=33.775867, alt=0, horizon=30)
        assert aos is None or aos > utc_time, "AOS may not occur with high horizon angle"


@pytest.mark.parametrize(
    ("utc_time", "lon", "lat", "alt", "horizon", "is_aos", "expect_result"),
    [
        # Basic pass over Atlanta
        (dt.datetime(2003, 3, 23, 0, 0, 0), -84.39733, 33.775867, 0, 0, True, True),
        # High elevation requirement (still detectable)
        (dt.datetime(2003, 3, 23, 0, 0, 0), -84.39733, 33.775867, 0, 30, True, True),
        # Near polar region (unlikely pass for ISS)
        (dt.datetime(2003, 3, 23, 0, 0, 0), 0, 89.0, 0, 0, True, False),
        # AOL detection instead of AOS
        (dt.datetime(2003, 3, 23, 0, 0, 0), -84.39733, 33.775867, 0, 0, False, True),
        # Satellite below horizon at start time (should still detect future AOS)
        (dt.datetime(2003, 3, 23, 0, 0, 0), -84.39733, 33.775867, 0, 10, True, True),
        # Observer at sea level vs. high altitude (compare AOS detectability)
        (dt.datetime(2003, 3, 23, 0, 0, 0), -84.39733, 33.775867, 3000, 10, True, False),
        # Horizon angle just below satellite max elevation (should fail)
        (dt.datetime(2003, 3, 23, 0, 0, 0), -84.39733, 33.775867, 0, 89, True, False),
        # Negative altitude (below sea level, e.g., Dead Sea)
        (dt.datetime(2003, 3, 23, 0, 0, 0), 35.5, 31.5, -430, 10, True, True),
        # Equator, zero horizon — should always detect AOS if satellite is visible
        (dt.datetime(2003, 3, 23, 0, 0, 0), 0, 0, 0, 0, True, True),
    ],
)
def test_crossing_detection(utc_time, lon, lat, alt, horizon, is_aos, expect_result):
    """Test horizon crossing detection under various observer conditions."""
    sat = Orbital(
        "ISS (ZARYA)",
        line1="1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
        line2="2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029",
    )

    result = (
        sat.find_aos(utc_time, lon, lat, alt, horizon) if is_aos else sat.find_aol(utc_time, lon, lat, alt, horizon)
    )

    if expect_result:
        assert result is not None, "Expected crossing, got None"
        assert result > utc_time, "Crossing must be after input time"
    else:
        assert result is None, "Expected no crossing, but got one"


@pytest.mark.parametrize(
    ("utc_time", "lon", "lat", "alt", "horizon", "is_aos", "expect_result"),
    [
        # Extreme case: ISS cannot reach 89° latitude with 80° elevation mask
        (dt.datetime(2003, 3, 23, 0, 0, 0), 0, 89.0, 0, 80, True, False),
    ],
)
def test_extreme_no_aos(utc_time, lon, lat, alt, horizon, is_aos, expect_result):
    """Test that no AOS is detected under extreme elevation and location constraints."""
    sat = Orbital(
        "ISS (ZARYA)",
        line1="1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
        line2="2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029",
    )

    result = (
        sat.find_aos(utc_time, lon, lat, alt, horizon) if is_aos else sat.find_aol(utc_time, lon, lat, alt, horizon)
    )

    if expect_result:
        assert result is not None, "Expected crossing, got None"
        assert result > utc_time, "Crossing must be after input time"
    else:
        assert result is None, "Expected no crossing, but got one"


@pytest.mark.parametrize(
    ("tle", "test_time", "expected_sign"),
    [
        # ISS: should be near equator, Z < 0 before ascending node
        (
            (
                "1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
                "2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029",
            ),
            dt.datetime(2003, 3, 23, 0, 0, 0),
            1,  # Expect Z velocity > 0 at ascending node
        ),
        # METOP-A: polar orbit, ascending node should be near equator
        (
            (
                "1 29499U 06044A   13060.48822809  .00000017  00000-0  27793-4 0  9819",
                "2 29499  98.6639 121.6164 0001449  71.9056  43.3132 14.21510544330271",
            ),
            dt.datetime(2013, 3, 1, 12, 0, 0),
            1,
        ),
        # SUOMI NPP: test near epoch
        (
            (
                "1 37849U 11061A   13061.24611272  .00000048  00000-0  43679-4 0  4334",
                "2 37849  98.7444   1.0588 0001264  63.8791 102.8546 14.19528338 69643",
            ),
            dt.datetime(2013, 3, 2, 6, 0, 0),
            1,
        ),
    ],
)
def test_get_last_an_time_velocity_sign(tle, test_time, expected_sign):
    """Test that get_last_an_time returns a time with positive Z velocity (ascending node)."""
    sat = Orbital("TestSat", line1=tle[0], line2=tle[1])
    an_time = sat.get_last_an_time(test_time)
    _, vel = sat.get_position(an_time, normalize=False)
    assert vel[2] * expected_sign > 0, f"Expected Z velocity sign {expected_sign}, got {vel[2]}"


@pytest.mark.parametrize(
    ("tle", "test_time"),
    [
        # ISS
        (
            (
                "1 25544U 98067A   03097.78853147  .00021906  00000-0  28403-3 0  8652",
                "2 25544  51.6361  13.7980 0004256  35.6671  59.2566 15.58778559250029",
            ),
            dt.datetime(2003, 3, 23, 0, 0, 0),
        ),
        # METOP-A
        (
            (
                "1 29499U 06044A   13060.48822809  .00000017  00000-0  27793-4 0  9819",
                "2 29499  98.6639 121.6164 0001449  71.9056  43.3132 14.21510544330271",
            ),
            dt.datetime(2013, 3, 1, 12, 0, 0),
        ),
    ],
)
def test_get_last_an_time_z_velocity_positive(tle, test_time):
    """Verify that Z velocity is positive at the computed ascending node time."""
    sat = Orbital("TestSat", line1=tle[0], line2=tle[1])
    an_time = sat.get_last_an_time(test_time)
    _, vel = sat.get_position(an_time, normalize=False)
    assert vel[2] > 0, f"Expected positive Z velocity, got {vel[2]}"


@pytest.mark.parametrize(
    ("tle", "test_time"),
    [
        # SUOMI NPP
        (
            (
                "1 37849U 11061A   13061.24611272  .00000048  00000-0  43679-4 0  4334",
                "2 37849  98.7444   1.0588 0001264  63.8791 102.8546 14.19528338 69643",
            ),
            dt.datetime(2013, 3, 2, 6, 0, 0),
        ),
    ],
)
def test_get_last_an_time_period_consistency(tle, test_time):
    """Ensure the time between consecutive ascending nodes matches the orbital period."""
    sat = Orbital("TestSat", line1=tle[0], line2=tle[1])
    an1 = sat.get_last_an_time(test_time)
    an2 = sat.get_last_an_time(an1 - np.timedelta64(10, "m"))
    period = (an1 - an2).astype("timedelta64[s]").item().total_seconds() / 60.0
    expected_period = 1440.0 / sat.tle.mean_motion  # minutes
    assert abs(period - expected_period) < 1.0, f"Period mismatch: got {period}, expected {expected_period}"


@pytest.mark.parametrize(
    ("tle", "test_time"),
    [
        # METOP-A near node
        (
            (
                "1 29499U 06044A   13060.48822809  .00000017  00000-0  27793-4 0  9819",
                "2 29499  98.6639 121.6164 0001449  71.9056  43.3132 14.21510544330271",
            ),
            dt.datetime(2013, 3, 1, 5, 30, 0),
        ),
    ],
)
def test_get_last_an_time_near_node(tle, test_time):
    """Confirm ascending node time is reasonably close to the input time."""
    sat = Orbital("TestSat", line1=tle[0], line2=tle[1])
    an_time = sat.get_last_an_time(test_time)
    delta = abs((np.datetime64(test_time) - an_time).astype("timedelta64[s]").item().total_seconds())
    assert delta < 7200, f"AN time too far from test time: {delta} seconds"


class TestGetObserverLook(unittest.TestCase):
    """Test the get_observer_look function."""

    def setUp(self):
        """Set up the test environment."""
        self.t = dt.datetime(2018, 1, 1, 0, 0, 0)
        self.sat_lon = np.array([[-89.5, -89.4, -89.5, -89.4], [-89.3, -89.2, -89.3, -89.2]])
        self.sat_lat = np.array([[45.5, 45.4, 45.5, 45.4], [45.3, 40.2, 45.3, 40.2]])
        self.sat_alt = 35786 * np.ones((2, 4))
        self.lon = np.array([[-85.5, -85.4, -89.5, -99.4], [-85.3, -89.2, -89.3, -79.2]])
        self.lat = np.array([[40.5, 40.4, 65.5, 45.4], [40.3, 40.2, 25.3, 40.2]])
        self.alt = np.zeros((2, 4))
        self.exp_azi = np.array(
            [
                [331.00275902, 330.95954165, 180, 86.435411],
                [330.91642994, 180, 0, 273.232073],
            ]
        )
        self.exp_elev = np.array(
            [
                [83.18070976, 83.17788976, 66.548467, 81.735221],
                [83.17507167, 90, 66.559906, 81.010018],
            ]
        )

    def test_basic_numpy(self):
        """Test with numpy array inputs."""
        azi, elev = get_observer_look(
            self.sat_lon,
            self.sat_lat,
            self.sat_alt,
            self.t,
            self.lon,
            self.lat,
            self.alt,
        )
        np.testing.assert_allclose(azi, self.exp_azi)
        np.testing.assert_allclose(elev, self.exp_elev)

    def test_basic_dask(self):
        """Test with dask array inputs."""
        sat_lon = da.from_array(self.sat_lon, chunks=2)
        sat_lat = da.from_array(self.sat_lat, chunks=2)
        sat_alt = da.from_array(self.sat_alt, chunks=2)
        lon = da.from_array(self.lon, chunks=2)
        lat = da.from_array(self.lat, chunks=2)
        alt = da.from_array(self.alt, chunks=2)
        azi, elev = get_observer_look(sat_lon, sat_lat, sat_alt, self.t, lon, lat, alt)
        np.testing.assert_allclose(azi.compute(), self.exp_azi)
        np.testing.assert_allclose(elev.compute(), self.exp_elev)

    def test_xarray_with_numpy(self):
        """Test with xarray DataArray with numpy array as inputs."""

        def _xarr_conv(input_array):
            return xr.DataArray(input_array)

        sat_lon = _xarr_conv(self.sat_lon)
        sat_lat = _xarr_conv(self.sat_lat)
        sat_alt = _xarr_conv(self.sat_alt)
        lon = _xarr_conv(self.lon)
        lat = _xarr_conv(self.lat)
        alt = _xarr_conv(self.alt)
        azi, elev = get_observer_look(sat_lon, sat_lat, sat_alt, self.t, lon, lat, alt)
        np.testing.assert_allclose(azi.data, self.exp_azi)
        np.testing.assert_allclose(elev.data, self.exp_elev)

    def test_xarray_with_dask(self):
        """Test with xarray DataArray with dask array as inputs."""

        def _xarr_conv(input_array):
            return xr.DataArray(da.from_array(input_array, chunks=2))

        sat_lon = _xarr_conv(self.sat_lon)
        sat_lat = _xarr_conv(self.sat_lat)
        sat_alt = _xarr_conv(self.sat_alt)
        lon = _xarr_conv(self.lon)
        lat = _xarr_conv(self.lat)
        alt = _xarr_conv(self.alt)
        azi, elev = get_observer_look(sat_lon, sat_lat, sat_alt, self.t, lon, lat, alt)
        np.testing.assert_allclose(azi.data.compute(), self.exp_azi)
        np.testing.assert_allclose(elev, self.exp_elev)

    def test_scalar(self):
        """Test with scalar inputs."""
        (azi, elev) = get_observer_look(0, 0, 30_000_000, self.t, 0, 0, 0)
        np.testing.assert_allclose(elev, 90)


class TestGetObserverLookNadir(unittest.TestCase):
    """Test the get_observer_look function when satellite is at nadir."""

    def setUp(self):
        """Set up for test observer at nadir.

        Note that rounding error differs between array types.
        With 1000 elements a test gives:
          1 error for basic numpy
          41 errors for basic dask
          63 errors for xarray with dask
          2 error for xarray with numpy
        """
        rng = np.random.RandomState(125)
        self.t = dt.datetime(2018, 1, 1, 0, 0, 0)
        self.sat_lon = 360 * rng.rand(100) - 180
        self.sat_lat = 180 * rng.rand(100) - 90
        self.sat_alt = rng.rand(100) + 850
        self.lon = self.sat_lon  # + 10E-17
        self.lat = self.sat_lat  # + 10E-17
        self.alt = np.zeros((100))
        self.exp_elev = np.zeros((100)) + 90

    def test_basic_numpy(self):
        """Test with numpy array inputs."""
        azi, elev = get_observer_look(
            self.sat_lon,
            self.sat_lat,
            self.sat_alt,
            self.t,
            self.lon,
            self.lat,
            self.alt,
        )
        assert np.sum(np.isnan(azi)) == 0
        assert not np.isnan(azi).any()
        np.testing.assert_allclose(elev, self.exp_elev)

    def test_basic_dask(self):
        """Test with dask array inputs."""
        sat_lon = da.from_array(self.sat_lon, chunks=2)
        sat_lat = da.from_array(self.sat_lat, chunks=2)
        sat_alt = da.from_array(self.sat_alt, chunks=2)
        lon = da.from_array(self.lon, chunks=2)
        lat = da.from_array(self.lat, chunks=2)
        alt = da.from_array(self.alt, chunks=2)
        azi, elev = get_observer_look(sat_lon, sat_lat, sat_alt, self.t, lon, lat, alt)
        assert np.sum(np.isnan(azi)) == 0
        assert not np.isnan(azi).any()
        np.testing.assert_allclose(elev.compute(), self.exp_elev)

    def test_xarray_with_numpy(self):
        """Test with xarray DataArray with numpy array as inputs."""

        def _xarr_conv(input_array):
            return xr.DataArray(input_array)

        sat_lon = _xarr_conv(self.sat_lon)
        sat_lat = _xarr_conv(self.sat_lat)
        sat_alt = _xarr_conv(self.sat_alt)
        lon = _xarr_conv(self.lon)
        lat = _xarr_conv(self.lat)
        alt = _xarr_conv(self.alt)
        azi, elev = get_observer_look(sat_lon, sat_lat, sat_alt, self.t, lon, lat, alt)
        assert np.sum(np.isnan(azi)) == 0
        assert not np.isnan(azi).any()
        np.testing.assert_allclose(elev.data, self.exp_elev)

    def test_xarray_with_dask(self):
        """Test with xarray DataArray with dask array as inputs."""

        def _xarr_conv(input_array):
            return xr.DataArray(da.from_array(input_array, chunks=2))

        sat_lon = _xarr_conv(self.sat_lon)
        sat_lat = _xarr_conv(self.sat_lat)
        sat_alt = _xarr_conv(self.sat_alt)
        lon = _xarr_conv(self.lon)
        lat = _xarr_conv(self.lat)
        alt = _xarr_conv(self.alt)
        azi, elev = get_observer_look(sat_lon, sat_lat, sat_alt, self.t, lon, lat, alt)
        assert np.sum(np.isnan(azi)) == 0
        assert not np.isnan(azi).any()
        np.testing.assert_allclose(elev, self.exp_elev)


class TestRegressions(unittest.TestCase):
    """Test regressions."""

    def test_63(self):
        """Check that no runtimewarning is raised, #63."""
        warnings.filterwarnings("error")
        orb = Orbital(
            "Suomi-NPP",
            line1="1 37849U 11061A   19292.84582509  .00000011  00000-0  25668-4 0  9997",
            line2="2 37849  98.7092 229.3263 0000715  98.5313 290.6262 14.19554485413345",
        )
        orb.get_next_passes(dt.datetime(2019, 10, 21, 16, 0, 0), 12, 123.29736, -13.93763, 0)
        warnings.filterwarnings("default")


@pytest.mark.parametrize(
    "dtime",
    [
        dt.datetime(2024, 6, 25, 11, 0, 18),
        dt.datetime(2024, 6, 25, 11, 5, 0, 0, dt.timezone.utc),
        np.datetime64("2024-06-25T11:10:00.000000"),
    ],
)
def test_get_last_an_time_scalar_input(dtime):
    """Test getting the time of the last ascending node - input time is a scalar."""
    orb = Orbital(
        "NOAA-20",
        line1="1 43013U 17073A   24176.73674251  .00000000  00000+0  11066-3 0 00014",
        line2="2 43013  98.7060 114.5340 0001454 139.3958 190.7541 14.19599847341971",
    )

    expected = np.datetime64("2024-06-25T10:44:18.234375")
    result = orb.get_last_an_time(dtime)
    assert abs(expected - result) < np.timedelta64(1, "s")


@pytest.mark.parametrize(
    "dtime",
    [
        dt.datetime(2024, 6, 25, 11, 5, 0, 0, dt.timezone(dt.timedelta(hours=1))),
    ],
)
def test_get_last_an_time_wrong_input(dtime):
    """Test getting the time of the last ascending node - wrong input."""
    orb = Orbital(
        "NOAA-20",
        line1="1 43013U 17073A   24176.73674251  .00000000  00000+0  11066-3 0 00014",
        line2="2 43013  98.7060 114.5340 0001454 139.3958 190.7541 14.19599847341971",
    )

    expected = "UTC time expected! Parsing a timezone aware datetime object requires it to be UTC!"
    with pytest.raises(ValueError, match=expected):
        _ = orb.get_last_an_time(dtime)
