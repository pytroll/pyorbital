#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014, 2022, 2024 Pytroll Community

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

"""Unit testing the Astronomy methods and functions."""


from datetime import datetime

import dask.array as da
import numpy as np
import numpy.typing as npt
import pytest

import pyorbital.astronomy as astr

try:
    from xarray import DataArray
except ImportError:
    DataArray = None


def _create_dask_array(input_list: list, dtype: npt.DTypeLike) -> da.Array:
    """Create a dummy dask array for testing."""
    np_arr = np.array(input_list, dtype=dtype)
    return da.from_array(np_arr)


def _create_xarray_numpy(input_list: list, dtype: npt.DTypeLike) -> DataArray:
    """Create a dummy xarray DataArray for testing."""
    np_arr = np.array(input_list, dtype=dtype)
    return DataArray(np_arr)


def _create_xarray_dask(input_list: list, dtype: npt.DTypeLike) -> DataArray:
    """Create a dummy daskified xarray DataArray for testing."""
    dask_arr = _create_dask_array(input_list, dtype)
    return DataArray(dask_arr)


class TestAstronomy:
    """Testing the Astronomy class."""

    @pytest.mark.parametrize(
        ("dt", "exp_jdays", "exp_j2000"),
        [
            (datetime(2000, 1, 1, 12, 0), 2451545.0, 0),
            (datetime(2009, 10, 8, 14, 30), 2455113.1041666665, 3568.1041666666665),
        ]
    )
    def test_jdays(self, dt, exp_jdays, exp_j2000):
        """Test julian day functions."""
        assert astr.jdays(dt) == exp_jdays
        assert astr.jdays2000(dt) == exp_j2000

    @pytest.mark.parametrize(
        ("lon", "lat", "exp_theta"),
        [
            # Norrkoping
            (16.1833, 58.6167, 60.371433482557833),
            (0.0, 0.0, 1.8751916863323426),
        ]
    )
    @pytest.mark.parametrize(
        ("dtype", "array_construct"),
        [
            (None, None),
            (np.float32, np.array),
            (np.float64, np.array),
            (np.float32, _create_dask_array),
            (np.float64, _create_dask_array),
            (np.float32, _create_xarray_numpy),
            (np.float64, _create_xarray_numpy),
            (np.float32, _create_xarray_dask),
            (np.float64, _create_xarray_dask),
        ]
    )
    def test_sunangles(self, lon, lat, exp_theta, dtype, array_construct):
        """Test the sun-angle calculations."""
        if array_construct is None and dtype is not None:
            pytest.skip(reason="Xarray dependency unavailable")

        time_slot = datetime(2011, 9, 23, 12, 0)
        abs_tolerance = 1e-8
        if dtype is not None:
            lon = array_construct([lon], dtype=dtype)
            lat = array_construct([lat], dtype=dtype)
            if np.dtype(dtype).itemsize < 8:
                abs_tolerance = 1e-4

        sun_theta = astr.sun_zenith_angle(time_slot, lon, lat)
        if dtype is None:
            assert sun_theta == pytest.approx(exp_theta, abs=abs_tolerance)
            assert isinstance(sun_theta, float)
        else:
            assert sun_theta.dtype == dtype
            np.testing.assert_allclose(sun_theta, exp_theta, atol=abs_tolerance)
            assert isinstance(sun_theta, type(lon))

    def test_sun_earth_distance_correction(self):
        """Test the sun-earth distance correction."""
        utc_time = datetime(2022, 6, 15, 12, 0, 0)
        corr = astr.sun_earth_distance_correction(utc_time)
        corr_exp = 1.0156952156742332
        assert corr == pytest.approx(corr_exp, abs=1e-8)
