"""Unit testing the Astronomy methods and functions."""


import datetime as dt

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
        ("dat", "exp_jdays", "exp_j2000"),
        [
            (dt.datetime(2000, 1, 1, 12, 0), 2451545.0, 0),
            (dt.datetime(2009, 10, 8, 14, 30), 2455113.1041666665, 3568.1041666666665),
        ]
    )
    def test_jdays(self, dat, exp_jdays, exp_j2000):
        """Test julian day functions."""
        assert astr.jdays(dat) == exp_jdays
        assert astr.jdays2000(dat) == exp_j2000

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

        time_slot = dt.datetime(2011, 9, 23, 12, 0)
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

    def test_sun_azimuth_angle(self):
        """Solar azimuth is exposed clockwise from north for array coordinates."""
        time_slot = dt.datetime(2020, 3, 20, 9, 0)
        azimuth = astr.sun_azimuth_angle(time_slot, np.array([0.0, 10.0]), np.array([0.0, 45.0]))

        np.testing.assert_allclose(azimuth, [89.87886063, 133.28324397], atol=1e-8)

    @pytest.mark.parametrize(
        ("lon", "lat"),
        [
            (10, 50),
            (0, 0),
            (16, 59),
            (-75, -30),
        ]
    )
    def test_sun_zenith_angle_integer_scalar(self, lon, lat):
        """Test that an integer scalar lon/lat neither crashes nor truncates to whole degrees."""
        # The #156 dtype-preservation recast read ``.dtype`` off the raw ``lon`` argument,
        # so an ``int`` scalar raised AttributeError instead of returning the angle.
        time_slot = dt.datetime(2011, 9, 23, 12, 0)
        expected = astr.sun_zenith_angle(time_slot, float(lon), float(lat))
        result = astr.sun_zenith_angle(time_slot, lon, lat)
        assert result == pytest.approx(expected, abs=1e-8)
        assert isinstance(result, float)

    @pytest.mark.parametrize("dtype", [np.int32, np.int64])
    def test_sun_zenith_angle_integer_array(self, dtype):
        """Test that integer-dtype arrays match the float result instead of truncating to whole degrees."""
        time_slot = dt.datetime(2011, 9, 23, 12, 0)
        lon = np.array([10, 20, 30, -75], dtype=dtype)
        lat = np.array([50, 55, 60, -30], dtype=dtype)
        result = astr.sun_zenith_angle(time_slot, lon, lat)
        expected = astr.sun_zenith_angle(time_slot, lon.astype(np.float64), lat.astype(np.float64))
        assert np.issubdtype(result.dtype, np.floating)
        np.testing.assert_allclose(result, expected, atol=1e-8)

    def test_sun_earth_distance_correction(self):
        """Test the sun-earth distance correction."""
        utc_time = dt.datetime(2022, 6, 15, 12, 0, 0)
        corr = astr.sun_earth_distance_correction(utc_time)
        corr_exp = 1.0156952156742332
        assert corr == pytest.approx(corr_exp, abs=1e-8)
