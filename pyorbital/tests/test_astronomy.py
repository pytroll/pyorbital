#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014, 2022 Pytroll Community

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

from datetime import datetime

import pytest

import pyorbital.astronomy as astr


class TestAstronomy:

    def test_jdays(self):
        """Test julian day functions."""
        t = datetime(2000, 1, 1, 12, 0)
        assert astr.jdays(t) == 2451545.0
        assert astr.jdays2000(t) == 0
        t = datetime(2009, 10, 8, 14, 30)
        assert astr.jdays(t) == 2455113.1041666665
        assert astr.jdays2000(t) == 3568.1041666666665

    def test_sunangles(self):
        """Test the sun-angle calculations."""
        lat, lon = 58.6167, 16.1833  # Norrkoping
        time_slot = datetime(2011, 9, 23, 12, 0)

        sun_theta = astr.sun_zenith_angle(time_slot, lon, lat)
        assert sun_theta == pytest.approx(60.371433482557833, abs=1e-8)
        sun_theta = astr.sun_zenith_angle(time_slot, 0., 0.)
        assert sun_theta == pytest.approx(1.8751916863323426, abs=1e-8)

    def test_sun_earth_distance_correction(self):
        """Test the sun-earth distance correction."""
        utc_time = datetime(2022, 6, 15, 12, 0, 0)
        corr = astr.sun_earth_distance_correction(utc_time)
        corr_exp = 1.0156952156742332
        assert corr == pytest.approx(corr_exp, abs=1e-8)
