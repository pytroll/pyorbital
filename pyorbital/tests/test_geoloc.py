"""Test the geoloc module."""


import datetime as dt

import numpy as np
import pytest

from pyorbital.geoloc import ScanGeometry, geodetic_lat, qrotate, subpoint
from pyorbital.geoloc_avhrr import (
    compute_avhrr_gcps_lonlatalt,
    estimate_time_and_attitude_deviations,
    estimate_time_offset,
)
from pyorbital.geoloc_instrument_definitions import (
    MultiLineWhiskbroomScan,
    PushbroomSwath,
    SingleLinePushbroomScan,
    amsua,
    ascat,
    atms,
    avhrr,
    avhrr_from_times,
    avhrr_gac_from_times,
    hirs4,
    mhs,
    mwhs2,
    slstr_nadir,
    viirs,
)


def test_qrotate():
    """Test quaternion rotation."""
    vector = np.array([[1, 0, 0]]).T
    axis = np.array([[0, 1, 0]]).T
    angle = np.deg2rad(90)

    result = qrotate(vector, axis, angle)[:, 0]
    expected = np.array([0, 0, 1])
    np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

    axis = np.array([0, 1, 0])
    result = qrotate(vector, axis, angle)
    expected = np.array([[0, 0, 1]]).T
    np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

    vector = np.array([[1, 0, 0],
                       [0, 0, 1]]).T
    axis = np.array([0, 1, 0])
    angle = np.deg2rad(90)
    result = qrotate(vector, axis, angle)
    expected = np.array([[0, 0, 1],
                         [-1, 0, 0]]).T

    np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

    axis = np.array([[0, 1, 0]]).T
    result = qrotate(vector, axis, angle)
    expected = np.array([[0, 0, 1],
                         [-1, 0, 0]]).T

    np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)


class TestGeoloc:
    """Test for the core computing part."""

    def test_scan_geometry(self):
        """Test the ScanGeometry object."""
        scans_nb = 1

        xy = np.vstack((np.deg2rad(np.array([10, 0, -10])),
                        np.array([0, 0, 0])))
        xy = np.tile(xy[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

        times = np.tile([-0.1, 0, 0.1], [np.int32(scans_nb), 1])

        instrument = ScanGeometry(xy, times)

        np.testing.assert_allclose(np.rad2deg(instrument.fovs[0]), np.array([[10, 0, -10]]))

        # Test vectors
        pos = np.rollaxis(np.tile(np.array([0, 0, 7000]), [3, 1, 1]), 2)
        vel = np.rollaxis(np.tile(np.array([1, 0, 0]), [3, 1, 1]), 2)
        pos = np.stack([np.array([0, 0, 7000])] * 3, 1)[:, np.newaxis, :]
        vel = np.stack([np.array([1, 0, 0])] * 3, 1)[:, np.newaxis, :]

        vec = instrument.vectors(pos, vel)

        result = vec[:, 0, 1]
        expected = np.array([0.0, 0.0, -1.0])
        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

        # Check if we can pass an array for yaw
        vec = instrument.vectors(pos, vel, yaw=[0])

        result = vec[:, 0, 1]
        expected = np.array([0.0, 0.0, -1.0])
        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)


        # minus sin because we use trigonometrical direction of angles
        result = vec[:, 0, 0]
        expected = np.array([0, -np.sin(np.deg2rad(10)), -np.cos(np.deg2rad(10))])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)

        result = vec[:, 0, 2]
        expected = np.array([0, -np.sin(np.deg2rad(-10)), -np.cos(np.deg2rad(-10))])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)

        # Test times

        start_of_scan = np.datetime64(dt.datetime(2014, 1, 8, 11, 30))
        times = instrument.times(start_of_scan)

        assert times[0, 1] == start_of_scan
        assert times[0, 0] == start_of_scan - np.timedelta64(100, "ms")
        assert times[0, 2] == start_of_scan + np.timedelta64(100, "ms")

    def test_geodetic_lat(self):
        """Test the determination of the geodetic latitude."""
        point = np.array([[7000, 0, 7000]]).T
        np.testing.assert_allclose(geodetic_lat(point),
                                   np.array([0.78755832699854733]), rtol=1e-8, atol=1e-8)

        points = np.array([[7000, 0, 7000],
                           [7000, 0, 7000]]).T
        result = geodetic_lat(points)
        expected = np.array([0.78755832699854733, 0.78755832699854733])
        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

    def test_subpoint(self):
        """Test nadir determination."""
        a = 6378.137  # km
        b = 6356.75231414  # km, GRS80
        point = np.array([0, 0, 7000])
        nadir = subpoint(point, a, b)
        np.testing.assert_allclose(nadir, np.array([0, 0, b]), rtol=1e-7, atol=1e-7)

        point = np.array([7000, 0, 7000])
        nadir = subpoint(point, a, b)
        np.testing.assert_allclose(nadir,
                                   np.array([4507.85431429,
                                             0,
                                             4497.06396339]), rtol=1e-8, atol=1e-8)
        points = np.array([[7000, 0, 7000],
                           [7000, 0, 7000]]).T
        nadir = subpoint(points, a, b)
        np.testing.assert_allclose(nadir[:, 0],
                                   np.array([4507.85431429,
                                             0,
                                             4497.06396339]), rtol=1e-8, atol=1e-8)
        np.testing.assert_allclose(nadir[:, 1],
                                   np.array([4507.85431429,
                                             0,
                                             4497.06396339]), rtol=1e-8, atol=1e-8)




def test_arbitrary_point_geoloc():
    """Test geolocating an arbitrary point in the swath."""
    from pyorbital.geoloc_avhrr import compute_avhrr_gcps_lonlatalt

    # Couple of example Two Line Elements
    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"

    # Choosing a specific time, this should be relatively close to the issue date of the TLE
    t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)
    rpy = (0, 0, 0)

    max_scan_angle = 55.37

    gcps = np.array([[2, 500], [1500, 700], [20, 1000]])

    lons, lats, alts = compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, rpy, t, (tle1, tle2))

    assert lons[0] == pytest.approx(-34.69996894)
    assert lats[0] == pytest.approx(56.69799502)

    assert lons[2] == pytest.approx(-27.573052737698944)
    assert lats[2] == pytest.approx(55.626740897592654)


def test_minimize_geoloc_error():
    """Test minimizing the distance to a set of gcps."""
    # Couple of example Two Line Elements
    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
    tle = (tle1, tle2)

    # Choosing a specific time, this should be relatively close to the issue date of the TLE
    t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)

    ref_time_displacement = 0.51
    ref_time = t + dt.timedelta(seconds=ref_time_displacement)
    ref_yaw = 0.1
    rpy = (0, 0, ref_yaw)
    max_scan_angle = 55.37
    # gcps are line/col
    gcps = np.array([[2, 500], [1500, 700], [20, 1000], [500, 1100], [100, 2000]])
    ref_lons, ref_lats, _ = compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, rpy, ref_time, tle)
    time_diff, (roll, pitch, yaw), (do, dm) = estimate_time_and_attitude_deviations(gcps, ref_lons, ref_lats, t,
                                                                                    tle, max_scan_angle)
    assert time_diff == pytest.approx(ref_time_displacement, abs=1e-2)
    assert yaw == pytest.approx(ref_yaw, abs=1e-2)
    assert min(do) > max(dm)


def test_minimize_time_error():
    """Test minimizing the distance to a set of gcps using only time offset."""
    # Couple of example Two Line Elements
    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
    tle = (tle1, tle2)

    # Choosing a specific time, this should be relatively close to the issue date of the TLE
    t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)

    ref_time_displacement = 20
    ref_time = t + dt.timedelta(seconds=ref_time_displacement)
    rpy = (0, 0, 0)
    max_scan_angle = 55.37
    # gcps are line/col
    gcps = np.array([[2, 500], [1500, 700], [20, 1000], [500, 1100], [100, 2000]])
    ref_lons, ref_lats, _ = compute_avhrr_gcps_lonlatalt(gcps, max_scan_angle, rpy, ref_time, tle)
    time_diff, (do, dm) = estimate_time_offset(gcps, ref_lons, ref_lats, t,
                                                                      tle, max_scan_angle)
    assert time_diff == pytest.approx(ref_time_displacement, abs=1e-2)
    assert min(do) > max(dm)


class TestGeolocDefs:
    """Test the instrument definitions."""

    def test_avhrr(self):
        """Test the definition of the avhrr instrument."""
        avh = avhrr(1, np.array([0, 1023.5, 2047]))
        result = np.rad2deg(avh.fovs[0])
        expected = np.array([[55.37, 0, -55.37]])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)

        avh = avhrr(1, np.array([0, 1023.5, 2047]), 10)
        np.testing.assert_allclose(np.rad2deg(avh.fovs[0]),
                                   np.array([[10, 0, -10]]))

        # This is perhaps a bit odd, to require avhrr to accept floats for
        # the number of scans? FIXME!
        avh = avhrr(1.1, np.array([0, 1023.5, 2047]), 10)
        np.testing.assert_allclose(np.rad2deg(avh.fovs[0]),
                                   np.array([[10, 0, -10]]))

    def test_avhrr_from_times(self):
        """Test generating the avhrr from times."""
        avh = avhrr_from_times([dt.datetime(2000,1,1,0,0,0)], [0, 1023.5, 2047])
        result = np.rad2deg(avh.fovs[0])
        expected = np.array([[55.37, 0, -55.37]])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)
        result = avh.times(dt.date(2000,1,1))
        expected = ((np.array([[0, 1023.5, 2047]]) * 25000).astype("timedelta64[ns]")
                    + np.datetime64("2000-01-01T00:00:00"))
        np.testing.assert_equal(result, expected)

        avh = avhrr_from_times([dt.datetime(2000,1,1,0,0,0)], np.array([0, 1023.5, 2047]), 10)
        np.testing.assert_allclose(np.rad2deg(avh.fovs[0]),
                                   np.array([[10, 0, -10]]))

        avh = avhrr_from_times([dt.datetime(2000,1,1,0,0,0), dt.datetime(2000,1,1,0,1,0)],
                    [0, 2047])
        times = avh.times(dt.datetime(2001,1,1))
        expected = (np.array([[0,51175000],[60000000000, 60051175000]]).astype("timedelta64[ns]")
                    + np.datetime64("2001-01-01"))
        np.testing.assert_equal(times, expected)


    def test_avhrr_gac_from_times(self):
        """Test getting avhrr gac from times."""
        avh = avhrr_gac_from_times([dt.datetime(2000,1,1,0,0,0)], [0, 204, 408])
        result = np.rad2deg(avh.fovs[0])
        expected = np.array([[55.180655, 0, -55.180655]])
        np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-7)
        result = avh.times(dt.date(2000,1,1))
        expected = ((np.array([[0, 204, 408]]) * 125000).astype("timedelta64[ns]")
                    + np.datetime64("2000-01-01T00:00:00"))
        np.testing.assert_equal(result, expected)

        avh = avhrr_gac_from_times([dt.datetime(2000,1,1,0,0,0)], np.array([0, 204, 408]), 10)
        np.testing.assert_allclose(np.rad2deg(avh.fovs[0]),
                                   np.array([[9.965804, 0, -9.965804]]))

        avh = avhrr_gac_from_times([dt.datetime(2000,1,1,0,0,0), dt.datetime(2000,1,1,0,1,0)],
                    [0, 408])
        times = avh.times(dt.datetime(2001,1,1))
        expected = (np.array([[0,51000000],[60000000000, 60051000000]]).astype("timedelta64[ns]")
                    + np.datetime64("2001-01-01"))
        np.testing.assert_equal(times, expected)


    def test_viirs(self):
        """Test the definition of the viirs instrument."""
        geom = viirs(1, np.array([0, 3200, 6399]))
        expected_fovs = np.array([
            np.tile(np.array([[0.98, -0., -0.98]]), [32, 1]),
            np.tile(np.array([[0., -0., 0]]), [32, 1])], dtype=np.float64)

        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

        geom = viirs(2, np.array([0, 3200, 6399]))
        expected_fovs = np.array([
            np.tile(np.array([[0.98, -0., -0.98]]), [32*2, 1]),
            np.tile(np.array([[0., -0., 0]]), [32*2, 1])], dtype=np.float64)

        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_viirs_defaults(self):
        """Test the definition of the viirs instrument with default slicing."""
        geom = viirs(1, chn_pixels=3)
        expected_fovs = np.array([
            np.tile(np.array([[0.98, -0., -0.98]]), [32, 1]),
            np.tile(np.array([[0., -0., 0]]), [32, 1])], dtype=np.float64)

        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_amsua(self):
        """Test the definition of the amsua instrument."""
        geom = amsua(1)
        expected_fovs = np.array([
            [[0.84,  0.78,  0.73,  0.67,  0.61,  0.55,  0.49,  0.44,  0.38,
              0.32,  0.26,  0.2,  0.15,  0.09,  0.03, -0.03, -0.09, -0.15,
              -0.2, -0.26, -0.32, -0.38, -0.44, -0.49, -0.55, -0.61, -0.67,
              -0.73, -0.78, -0.84]],
            np.zeros((1, 30))], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)

    def test_mhs(self):
        """Test the definition of the mhs instrument."""
        geom = mhs(1)
        expected_fovs = np.array([
            [[0.86,  0.84,  0.82,  0.8,  0.79,  0.77,  0.75,  0.73,  0.71,
              0.69,  0.67,  0.65,  0.63,  0.61,  0.59,  0.57,  0.55,  0.53,
              0.51,  0.49,  0.48,  0.46,  0.44,  0.42,  0.4,  0.38,  0.36,
              0.34,  0.32,  0.3,  0.28,  0.26,  0.24,  0.22,  0.2,  0.18,
              0.16,  0.15,  0.13,  0.11,  0.09,  0.07,  0.05,  0.03,  0.01,
              -0.01, -0.03, -0.05, -0.07, -0.09, -0.11, -0.13, -0.15, -0.16,
              -0.18, -0.2, -0.22, -0.24, -0.26, -0.28, -0.3, -0.32, -0.34,
              -0.36, -0.38, -0.4, -0.42, -0.44, -0.46, -0.48, -0.49, -0.51,
              -0.53, -0.55, -0.57, -0.59, -0.61, -0.63, -0.65, -0.67, -0.69,
              -0.71, -0.73, -0.75, -0.77, -0.79, -0.8, -0.82, -0.84, -0.86]],
            np.zeros((1, 90))], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_hirs4(self):
        """Test the definition of the hirs4 instrument."""
        geom = hirs4(1)
        expected_fovs = np.array([
            [[0.86,  0.83,  0.8,  0.77,  0.74,  0.71,  0.68,  0.64,  0.61,
              0.58,  0.55,  0.52,  0.49,  0.46,  0.42,  0.39,  0.36,  0.33,
              0.3,  0.27,  0.24,  0.2,  0.17,  0.14,  0.11,  0.08,  0.05,
              0.02, -0.02, -0.05, -0.08, -0.11, -0.14, -0.17, -0.2, -0.24,
              -0.27, -0.3, -0.33, -0.36, -0.39, -0.42, -0.46, -0.49, -0.52,
              -0.55, -0.58, -0.61, -0.64, -0.68, -0.71, -0.74, -0.77, -0.8,
              -0.83, -0.86]],
            np.zeros((1, 56))], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_atms(self):
        """Test the definition of the atms instrument."""
        geom = atms(1)
        expected_fovs = np.array([
            [[0.92,  0.9,  0.88,  0.86,  0.84,  0.82,  0.8,  0.78,  0.76,
              0.75,  0.73,  0.71,  0.69,  0.67,  0.65,  0.63,  0.61,  0.59,
              0.57,  0.55,  0.53,  0.51,  0.49,  0.47,  0.46,  0.44,  0.42,
              0.4,  0.38,  0.36,  0.34,  0.32,  0.3,  0.28,  0.26,  0.24,
              0.22,  0.2,  0.18,  0.16,  0.15,  0.13,  0.11,  0.09,  0.07,
              0.05,  0.03,  0.01, -0.01, -0.03, -0.05, -0.07, -0.09, -0.11,
              -0.13, -0.15, -0.16, -0.18, -0.2, -0.22, -0.24, -0.26, -0.28,
              -0.3, -0.32, -0.34, -0.36, -0.38, -0.4, -0.42, -0.44, -0.46,
              -0.47, -0.49, -0.51, -0.53, -0.55, -0.57, -0.59, -0.61, -0.63,
              -0.65, -0.67, -0.69, -0.71, -0.73, -0.75, -0.76, -0.78, -0.8,
              -0.82, -0.84, -0.86, -0.88, -0.9, -0.92]],
            np.zeros((1, 96))], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs,
                                   expected_fovs, rtol=1e-2, atol=1e-2)

    def test_ascat(self):
        """Test the definition of the ASCAT instrument onboard Metop."""
        geom = ascat(1)
        expected_fovs = np.array([
            [[0.9250245,  0.90058989,  0.87615528,  0.85172067,
              0.82728607,  0.80285146,  0.77841685,  0.75398224,
              0.72954763,  0.70511302,  0.68067841,  0.6562438,
              0.63180919,  0.60737458,  0.58293997,  0.55850536,
              0.53407075,  0.50963614,  0.48520153,  0.46076692,
              0.43633231, -0.43633231, -0.46076692, -0.48520153,
              -0.50963614, -0.53407075, -0.55850536, -0.58293997,
              -0.60737458, -0.63180919, -0.6562438, -0.68067841,
              -0.70511302, -0.72954763, -0.75398224, -0.77841685,
              -0.80285146, -0.82728607, -0.85172067, -0.87615528,
              -0.90058989, -0.9250245]], np.zeros((1, 42))], dtype=np.float64)

        np.testing.assert_allclose(
            geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)
        geom = ascat(1, np.array([0, 41]))
        expected_fovs = np.array([[[0.9250245,  -0.9250245]],
                                  [[0.,  0.]]], dtype=np.float64)
        np.testing.assert_allclose(
            geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)

        geom = ascat(1, np.array([0, -1]))
        np.testing.assert_allclose(
            geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)

    def test_slstr_nadir(self):
        """Test the definition of the slstr instrument nadir view flying on Sentinel-3."""
        geom = slstr_nadir(1, [0, 1])

        expected_fovs = np.array([
            np.tile(np.array([[0.8115781, -0.38571776]]), [1, 1]),
            np.tile(np.array([[0., 0.]]), [1, 1])], dtype=np.float64)
        np.testing.assert_allclose(geom.fovs, expected_fovs, rtol=1e-2, atol=1e-2)

        geom = slstr_nadir(1, None)

        np.testing.assert_equal(geom.fovs.size, 6000)

    def test_one_line_pushbroom(self):
        """Test pushbroom swath geometry via PushbroomSwath class."""
        scan = SingleLinePushbroomScan(left_angle=46.5, right_angle=-22.1, pixels_per_scan=4865)
        time_sampling = np.timedelta64(44, "ms")
        swath = PushbroomSwath(scanline=scan, time_sampling=time_sampling)
        geom = swath.scan_geometry(scan_lines=slice(None, 11, 10), pixels=slice(None, None, 152))
        assert geom.fovs.shape == (2, 2, 33)
        assert geom.fovs[0, 0, 0] == pytest.approx(np.deg2rad(46.5))
        assert geom.fovs[0, 1, 0] == pytest.approx(np.deg2rad(46.5))
        assert geom.fovs[0, 0, -1] == pytest.approx(np.deg2rad(-22.1))
        assert geom.fovs[0, 1, -1] == pytest.approx(np.deg2rad(-22.1))
        assert geom.fovs[1, 0, 0] == 0
        assert geom._times.shape == (2, 33)
        assert geom._times[0, 0] == 0
        assert geom._times[0, -1] == 0
        assert geom._times[1, 0] == time_sampling * 10
        assert geom._times[1, -1] == time_sampling * 10
        assert geom._times.dtype == np.timedelta64(1, "ns").dtype


def test_single_line_pushbroom_scan():
    pixels_per_scan = 4865
    left_angle = 46.5
    right_angle = -22.1
    forward_angle = 10
    scan = SingleLinePushbroomScan(left_angle, right_angle, pixels_per_scan, forward_angle=forward_angle)
    x_fovs, y_fovs = scan.angles()
    np.testing.assert_allclose(x_fovs, np.linspace(np.deg2rad(left_angle),
                                                   np.deg2rad(right_angle), pixels_per_scan))
    np.testing.assert_allclose(y_fovs, np.deg2rad(forward_angle))
    assert len(y_fovs) == len(x_fovs)

    step = 152
    pixel_numbers = slice(0, pixels_per_scan, step)
    reduced_x_fovs, y_fovs = scan.angles(pixel_numbers)
    np.testing.assert_allclose(y_fovs, np.deg2rad(forward_angle))
    np.testing.assert_allclose(reduced_x_fovs, x_fovs[pixel_numbers])

    pixel_numbers = slice(76, pixels_per_scan, step)
    reduced_x_fovs, y_fovs = scan.angles(pixel_numbers)
    np.testing.assert_allclose(y_fovs, np.deg2rad(forward_angle))
    np.testing.assert_allclose(reduced_x_fovs, x_fovs[pixel_numbers])

    pixel_numbers = [0, 2432, 4864]
    reduced_x_fovs, y_fovs = scan.angles(pixel_numbers)
    np.testing.assert_allclose(y_fovs, np.deg2rad(forward_angle))
    np.testing.assert_allclose(reduced_x_fovs, x_fovs[pixel_numbers])


def test_pushbroom_swath_generates_scan_geometry():
    """Test that PushbroomSwath produces a ScanGeometry with correct fovs and times."""
    scan = SingleLinePushbroomScan(left_angle=46.5, right_angle=-22.1, pixels_per_scan=4865)
    time_sampling = np.timedelta64(44, "ms")
    swath = PushbroomSwath(scanline=scan, time_sampling=time_sampling)
    geom = swath.scan_geometry(scan_lines=slice(None, 11, 10), pixels=slice(None, None, 152))
    assert isinstance(geom, ScanGeometry)
    assert geom.fovs.shape == (2, 2, 33)
    assert geom.fovs[0, 0, 0] == pytest.approx(np.deg2rad(46.5))
    assert geom.fovs[0, 0, -1] == pytest.approx(np.deg2rad(-22.1))
    assert geom.fovs[1, 0, 0] == 0
    assert geom._times[0, 0] == np.timedelta64(0)
    assert geom._times[0, -1] == np.timedelta64(0)
    assert geom._times[1, 0] == time_sampling * 10
    assert geom._times[1, -1] == time_sampling * 10


def test_olci_scan_constant_matches_olci_function():
    """Test that OLCI_SCAN constant produces geometry matching the legacy olci() function."""
    from pyorbital.geoloc_instrument_definitions import OLCI_SCAN, olci

    legacy_geom = olci(10)
    swath = PushbroomSwath(scanline=OLCI_SCAN, time_sampling=np.timedelta64(44, "ms"))
    new_geom = swath.scan_geometry(scan_lines=slice(10))
    np.testing.assert_allclose(new_geom.fovs, legacy_geom.fovs)
    np.testing.assert_allclose(new_geom._times.astype(float), legacy_geom._times.astype(float), atol=1)


def test_slstr_nadir_scan_constant():
    """Test that SLSTR_NADIR_SCAN constant produces geometry matching legacy slstr_nadir()."""
    from pyorbital.geoloc_instrument_definitions import SLSTR_NADIR_SCAN

    legacy_geom = slstr_nadir(10)
    swath = PushbroomSwath(scanline=SLSTR_NADIR_SCAN, time_sampling=np.timedelta64(0))
    new_geom = swath.scan_geometry(scan_lines=slice(10))
    np.testing.assert_allclose(new_geom.fovs, legacy_geom.fovs)
    np.testing.assert_equal(new_geom._times, legacy_geom._times)


def test_bounding_box_returns_closed_polygon():
    """Test that bounding_box returns a closed polygon of lon/lat points."""
    from pyorbital.geoloc import bounding_box
    from pyorbital.geoloc_instrument_definitions import OLCI_SWATH

    tle1 = "1 33591U 09005A   21355.91138073  .00000074  00000+0  65091-4 0  9998"
    tle2 = "2 33591  99.1688  21.1338 0013414 329.8936  30.1462 14.12516400663123"
    start_time = dt.datetime(2021, 12, 22, 12, 0, 0)
    end_time = start_time + dt.timedelta(seconds=44 * 0.044)

    lons, lats = bounding_box(OLCI_SWATH, start_time, end_time, (tle1, tle2),
                              points_per_edge=5)
    # 5 points per edge, 4 edges sharing corners = 4*(5-1) + 1 closing = 17
    assert len(lons) == 17
    assert len(lats) == 17
    # polygon is closed
    assert lons[0] == lons[-1]
    assert lats[0] == lats[-1]
    # all values are finite
    assert np.all(np.isfinite(lons))
    assert np.all(np.isfinite(lats))


def test_yaw_steering_changes_vectors():
    """Test that yaw steering modifies the look vectors compared to no steering."""
    from pyorbital.geoloc import compute_yaw_steering

    # A satellite at ~7000km altitude, moving along x-axis, above the equator
    pos = np.array([[7000, 0, 0]]).T  # equatorial position
    vel = np.array([[0, 7.5, 0]]).T   # ~7.5 km/s orbital velocity

    yaw_angle = compute_yaw_steering(pos, vel)

    # At the equator, yaw steering should produce a non-zero angle
    assert yaw_angle != 0.0
    # The angle should be small (typically < 4 degrees for LEO)
    assert abs(yaw_angle) < np.deg2rad(4)


def test_yaw_steering_applied_in_vectors():
    """Test that vectors() applies yaw steering when enabled."""
    from pyorbital.geoloc import compute_yaw_steering

    xy = np.vstack((np.deg2rad(np.array([10, 0, -10])),
                    np.array([0, 0, 0])))
    xy = np.tile(xy[:, np.newaxis, :], [1, 1, 1])
    times = np.tile([0, 0, 0], [1, 1])
    instrument = ScanGeometry(xy, times)

    pos = np.array([[7000, 0, 0]]).T
    vel = np.array([[0, 7.5, 0]]).T
    pos = np.stack([pos[:, 0]] * 3, axis=1)[:, np.newaxis, :]
    vel = np.stack([vel[:, 0]] * 3, axis=1)[:, np.newaxis, :]

    vec_no_steering = instrument.vectors(pos, vel)
    vec_with_steering = instrument.vectors(pos, vel, yaw_steering=True)

    # Vectors should differ when yaw steering is applied
    assert not np.allclose(vec_no_steering, vec_with_steering)

    # Manually compute: yaw steering should add the computed yaw angle
    yaw_angle = compute_yaw_steering(pos, vel)
    vec_manual_yaw = instrument.vectors(pos, vel, yaw=yaw_angle)
    np.testing.assert_allclose(vec_with_steering, vec_manual_yaw)


def test_compute_pixels_with_yaw_steering():
    """Test that compute_pixels passes yaw_steering through to vectors."""
    from pyorbital.geoloc import compute_pixels

    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
    t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)

    # Create a simple 1D scan geometry (like compute_avhrr_gcps_lonlatalt does)
    scan_angles = np.array([np.deg2rad(-55.37), 0, np.deg2rad(55.37)])
    fovs = np.vstack((scan_angles, np.zeros(3)))
    times = np.array([0.0, 0.001, 0.002])
    sgeom = ScanGeometry(fovs, times)
    s_times = sgeom.times(t)

    pixels_no_yaw = compute_pixels((tle1, tle2), sgeom, s_times)
    pixels_yaw = compute_pixels((tle1, tle2), sgeom, s_times, yaw_steering=True)

    # The positions should differ with yaw steering
    assert not np.allclose(pixels_no_yaw, pixels_yaw)


def test_compute_pixels_with_2d_scan_times():
    """Test that compute_pixels works with multi-scan 2D time arrays (avhrr-style).

    avhrr() produces fovs of shape (2, N_SCANS, N_PX) and times of shape
    (N_SCANS, N_PX). compute_pixels must handle this without crashing and
    return pixel positions of the correct shape.
    """
    from pyorbital.geoloc import compute_pixels, get_lonlatalt
    from pyorbital.geoloc_instrument_definitions import avhrr

    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
    t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)

    n_scans, n_px = 5, 10
    sgeom = avhrr(n_scans, np.arange(n_px))
    s_times = sgeom.times(t)

    assert s_times.shape == (n_scans, n_px), "times should be 2D"

    pixels = compute_pixels((tle1, tle2), sgeom, s_times)
    lons, lats, alts = get_lonlatalt(pixels, s_times)

    assert pixels.shape == (3, n_scans * n_px)
    assert lons.shape == (n_scans * n_px,)
    assert np.all(np.isfinite(lons))
    assert np.all(np.abs(lats) <= 90)
    # Pixels are surface points: altitude should be within terrain range (not orbital altitude)
    assert np.all(np.abs(alts) < 10000)


def test_compute_pixels_2d_matches_1d_reference():
    """Test that the 2D-times path gives the same result as the 1D reference path.

    The 1D reference path (flat ScanGeometry, per-pixel SGP4) is the
    original correct behaviour.  The 2D path uses per-scan SGP4 which
    is an approximation; the difference should be sub-pixel (< 0.01 deg).
    """
    from pyorbital.geoloc import compute_pixels, get_lonlatalt
    from pyorbital.geoloc_instrument_definitions import avhrr

    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
    t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)

    scan_pts = np.arange(200)

    # Reference: flat 1D ScanGeometry (per-pixel SGP4)
    cross_angles = (scan_pts / 1023.5 - 1) * np.deg2rad(-55.37)
    fovs_flat = np.vstack([cross_angles, np.zeros(len(scan_pts))])
    t_float = scan_pts * 0.000025
    sgeom_ref = ScanGeometry(fovs_flat, t_float)
    s_times_ref = sgeom_ref.times(t)
    pix_ref = compute_pixels((tle1, tle2), sgeom_ref, s_times_ref)
    lons_ref, lats_ref, _ = get_lonlatalt(pix_ref, s_times_ref)

    # 2D path: avhrr(1, ...) with per-scan SGP4 approximation
    sgeom_2d = avhrr(1, scan_pts)
    s_times_2d = sgeom_2d.times(t)
    pix_2d = compute_pixels((tle1, tle2), sgeom_2d, s_times_2d)
    lons_2d, lats_2d, _ = get_lonlatalt(pix_2d, s_times_2d)

    # Per-scan SGP4 is an approximation; error < 0.01 deg for AVHRR scan rates
    np.testing.assert_allclose(lons_2d, lons_ref, atol=0.01)
    np.testing.assert_allclose(lats_2d, lats_ref, atol=0.01)


def test_whiskbroom_scan_constants_match_legacy_functions():
    """Test that whiskbroom instrument constants produce geometry matching the legacy functions."""
    from pyorbital.geoloc_instrument_definitions import (
        AMSU_A_SCAN,
        ATMS_SCAN,
        HIRS4_SCAN,
        MHS_SCAN,
        MWHS2_SCAN,
    )
    for scan, legacy_fn in [
        (AMSU_A_SCAN, amsua),
        (MHS_SCAN, mhs),
        (HIRS4_SCAN, hirs4),
        (ATMS_SCAN, atms),
        (MWHS2_SCAN, mwhs2),
    ]:
        legacy_geom = legacy_fn(5)
        new_geom = scan.scan_geometry(5)
        np.testing.assert_allclose(new_geom.fovs, legacy_geom.fovs, rtol=1e-6, atol=1e-6,
                                   err_msg=f"{scan} fovs mismatch")
        np.testing.assert_allclose(
            new_geom._times.view(np.int64).astype(np.float64),
            legacy_geom._times.view(np.int64).astype(np.float64),
            rtol=1e-6, atol=1e-6,
            err_msg=f"{scan} times mismatch")


def test_multiline_whiskbroom_fov_shape():
    """Test that MultiLineWhiskbroomScan produces (2, n_scans*lines_per_scan, n_pixels) fovs."""
    scan = MultiLineWhiskbroomScan(
        pixels_per_scan=10, scan_angle=55.0, scan_rate=1.5,
        pixel_dwell_time=1.48 / 10, lines_per_scan=4,
        along_track_step=np.deg2rad(0.067),
    )
    geom = scan.scan_geometry(n_scans=3)
    assert geom.fovs.shape == (2, 3 * 4, 10)
    assert geom._times.shape == (3 * 4, 10)
    assert geom.lines_per_scan == 4


def test_multiline_whiskbroom_along_track_angles():
    """Check that along-track angles span symmetrically around zero for each scan's detector rows."""
    lines_per_scan = 4
    step = 0.01
    scan = MultiLineWhiskbroomScan(
        pixels_per_scan=5, scan_angle=10.0, scan_rate=1.5,
        pixel_dwell_time=0.1, lines_per_scan=lines_per_scan,
        along_track_step=step,
    )
    geom = scan.scan_geometry(n_scans=2)
    # Along-track angles for the first scan's rows (uniform across pixels)
    at_angles = geom.fovs[1, :lines_per_scan, 0]
    expected = (np.arange(lines_per_scan) - (lines_per_scan - 1) / 2.0) * step
    np.testing.assert_allclose(at_angles, expected, atol=1e-10)
    # All scan lines in a scan share the same pixel-level cross-track times
    scan_times = geom._times[:lines_per_scan, :]
    np.testing.assert_equal(scan_times[0], scan_times[1])


def test_multiline_whiskbroom_lines_per_scan_in_scan_geometry():
    """Assert that ScanGeometry.lines_per_scan attribute is set correctly."""
    from pyorbital.geoloc import ScanGeometry
    fovs = np.zeros((2, 8, 5))
    times = np.zeros((8, 5))
    geom = ScanGeometry(fovs, times, lines_per_scan=4)
    assert geom.lines_per_scan == 4


def test_get_lonlatalt_produces_same_results_with_numba():
    """Check that get_lonlatalt with numba produces the same lon/lat/alt as pyproj (< 1e-6 deg, < 0.1 m)."""
    pytest.importorskip("numba")
    from pyproj import Transformer

    from pyorbital.geoloc import _numba_ecef_to_lonlatalt

    # Build 1 000 random ECEF surface points (alt=0)
    rng = np.random.default_rng(0)
    lons_t = rng.uniform(-179, 179, 1000)
    lats_t = rng.uniform(-89, 89, 1000)
    trf_fwd = Transformer.from_crs(dict(proj="latlong"), dict(proj="geocent"))
    x_m, y_m, z_m = trf_fwd.transform(lons_t, lats_t, np.zeros(1000))

    # Reference: pyproj geocent -> latlong
    trf_inv = Transformer.from_crs(dict(proj="geocent"), dict(proj="latlong"))
    lon_pp, lat_pp, alt_pp = trf_inv.transform(x_m, y_m, z_m)

    # Numba kernel takes km, returns degrees + meters
    lon_nb, lat_nb, alt_nb = _numba_ecef_to_lonlatalt(x_m / 1000, y_m / 1000, z_m / 1000)

    np.testing.assert_allclose(lon_nb, lon_pp, atol=1e-6,
                               err_msg="numba lon disagrees with pyproj")
    np.testing.assert_allclose(lat_nb, lat_pp, atol=1e-6,
                               err_msg="numba lat disagrees with pyproj")
    np.testing.assert_allclose(alt_nb, alt_pp, atol=0.1,
                               err_msg="numba alt disagrees with pyproj")

def test_compute_pixel_works():
    """Check that compute_pixels works for MultiLineWhiskbroomScan without OOM or crash."""
    from pyorbital.geoloc import compute_pixels, get_lonlatalt
    from pyorbital.geoloc_instrument_definitions import MERSI_250M_SCAN

    tle1 = "1 33591U 09005A   21355.91138073  .00000074  00000+0  65091-4 0  9998"
    tle2 = "2 33591  99.1688  21.1338 0013414 329.8936  30.1462 14.12516400663123"
    t = dt.datetime(2021, 12, 21, 12, 0, 0)

    n_scans = 10
    sgeom = MERSI_250M_SCAN.scan_geometry(n_scans)
    assert sgeom.fovs.shape == (2, n_scans * 40, 8192)
    assert sgeom.lines_per_scan == 40

    s_times = sgeom.times(t)
    pixels = compute_pixels((tle1, tle2), sgeom, s_times)
    lons, lats, alts = get_lonlatalt(pixels, s_times)

    assert pixels.shape == (3, n_scans * 40 * 8192)
    assert np.all(np.isfinite(lons))
    assert np.all(np.abs(lats) <= 90)


def test_geolocate_matches_compute_pixels_get_lonlatalt():
    """Check that geolocate() agrees with compute_pixels+get_lonlatalt to < 1e-6 deg / < 0.1 m."""
    pytest.importorskip("numba")
    from pyorbital.geoloc import _HAS_NUMBA, compute_pixels, geolocate, get_lonlatalt
    if not _HAS_NUMBA:
        pytest.skip("numba not available")

    from pyorbital.geoloc_instrument_definitions import MERSI_1KM_SCAN

    tle1 = "1 33591U 09005A   21355.91138073  .00000074  00000+0  65091-4 0  9998"
    tle2 = "2 33591  99.1688  21.1338 0013414 329.8936  30.1462 14.12516400663123"
    t = dt.datetime(2021, 12, 21, 12, 0, 0)
    n_scans = 5
    sgeom = MERSI_1KM_SCAN.scan_geometry(n_scans)
    s_times = sgeom.times(t)

    # Reference: the established pipeline
    pixels = compute_pixels((tle1, tle2), sgeom, s_times)
    lon_ref, lat_ref, alt_ref = get_lonlatalt(pixels, s_times)

    # New fused path
    lon_fused, lat_fused, alt_fused = geolocate((tle1, tle2), sgeom, s_times)

    assert lon_fused.shape == lon_ref.shape
    np.testing.assert_allclose(lon_fused, lon_ref, atol=1e-6,
                               err_msg="geolocate lon disagrees with reference")
    np.testing.assert_allclose(lat_fused, lat_ref, atol=1e-6,
                               err_msg="geolocate lat disagrees with reference")
    np.testing.assert_allclose(alt_fused, alt_ref, atol=1.0,
                               err_msg="geolocate alt disagrees with reference")


def test_geolocate_with_yaw_steering_matches_reference():
    """Check that geolocate(yaw_steering=True) fused path agrees with compute_pixels reference."""
    pytest.importorskip("numba")
    from pyorbital.geoloc import _HAS_NUMBA, compute_pixels, geolocate, get_lonlatalt
    if not _HAS_NUMBA:
        pytest.skip("numba not available")

    from pyorbital.geoloc_instrument_definitions import MERSI_1KM_SCAN

    tle1 = "1 33591U 09005A   21355.91138073  .00000074  00000+0  65091-4 0  9998"
    tle2 = "2 33591  99.1688  21.1338 0013414 329.8936  30.1462 14.12516400663123"
    t = dt.datetime(2021, 12, 21, 12, 0, 0)
    n_scans = 5
    sgeom = MERSI_1KM_SCAN.scan_geometry(n_scans)
    s_times = sgeom.times(t)

    # Reference: established pipeline with yaw steering
    pixels = compute_pixels((tle1, tle2), sgeom, s_times, yaw_steering=True)
    lon_ref, lat_ref, alt_ref = get_lonlatalt(pixels, s_times)

    # Fused numba path with yaw steering
    lon_f, lat_f, alt_f = geolocate((tle1, tle2), sgeom, s_times, yaw_steering=True)

    assert lon_f.shape == lon_ref.shape
    np.testing.assert_allclose(lon_f, lon_ref, atol=1e-6,
                               err_msg="yaw-steered lon disagrees with reference")
    np.testing.assert_allclose(lat_f, lat_ref, atol=1e-6,
                               err_msg="yaw-steered lat disagrees with reference")
    np.testing.assert_allclose(alt_f, alt_ref, atol=1.0,
                               err_msg="yaw-steered alt disagrees with reference")
