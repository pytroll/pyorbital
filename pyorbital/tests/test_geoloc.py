"""Test the geoloc module."""


import contextlib
import datetime as dt
import warnings

import numpy as np
import pytest

from pyorbital.config import config
from pyorbital.geoloc import (
    ScanGeometry,
    compute_pixels,
    geodetic_lat,
    geolocate,
    qrotate,
    subpoint,
)
from pyorbital.geoloc_avhrr import (
    compute_avhrr_gcps_lonlatalt,
    estimate_time_and_attitude_deviations,
    estimate_time_offset,
)
from pyorbital.geoloc_instrument_definitions import (
    FocalPlaneWhiskbroomScan,
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
from pyorbital.orbital import Orbital

# For tests whose assertions are convention-invariant (shapes, finiteness, or two
# code paths agreeing): the checks must hold whichever convention is selected, and
# relying on the default must still warn.
NADIR_CONVENTION_CASES = [
    pytest.param({}, True, id="default-warns"),
    pytest.param({"nadir_convention": "legacy"}, False, id="explicit-legacy"),
    pytest.param({"nadir_convention": "geocentric"}, False, id="geocentric"),
]

NADIR_CONVENTIONS = [
    pytest.param({}, True, "legacy", id="default-warns"),
    pytest.param({"nadir_convention": "legacy"}, False, "legacy", id="explicit-legacy"),
    pytest.param({"nadir_convention": "geocentric"}, False, "geocentric", id="geocentric"),
]


# For tests that reach the geolocation through APIs which take no convention
# argument (geoloc_avhrr, bounding_box): the choice can only be made
# process-wide, which is precisely what the config exists for.
NADIR_CONFIG_CASES = [
    pytest.param(None, True, id="default-warns"),
    pytest.param("legacy", False, id="config-legacy"),
    pytest.param("geocentric", False, id="config-geocentric"),
]


def _configured_convention(convention):
    """Context setting the process-wide convention, or leaving it unset."""
    if convention is None:
        return contextlib.nullcontext()
    return config.set(nadir_convention=convention)


def _expect_convention_warning(warns):
    """Return a context asserting the legacy-convention deprecation fires, or does not.

    Matching on the message keeps the assertion specific: an unrelated
    DeprecationWarning from elsewhere must not satisfy these tests.
    """
    if warns:
        return pytest.warns(DeprecationWarning, match="legacy nadir convention")
    return contextlib.nullcontext()


def _inclined_orbit_state():
    """A satellite position and a matching inclined-orbit velocity.

    The velocity must not be coplanar with the meridian: the local frame
    re-orthogonalises the nadir against it, and in the coplanar case that
    projects every convention onto -pos/|pos|, making them indistinguishable.
    """
    pos = np.array([[4507.0], [0.0], [5000.0]])
    orbit_normal = np.array([0.0, -np.sin(np.deg2rad(98.7)), np.cos(np.deg2rad(98.7))])
    vel = np.cross(orbit_normal, pos[:, 0]).reshape(3, 1)
    return pos, vel / np.sqrt((vel**2).sum()) * 7.5



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

    @pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
    def test_scan_geometry(self, kwargs, warns):
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

        context = _expect_convention_warning(warns)
        with context:
            vec = instrument.vectors(pos, vel, **kwargs)

        result = vec[:, 0, 1]
        expected = np.array([0.0, 0.0, -1.0])
        np.testing.assert_allclose(result, expected, rtol=1e-8, atol=1e-8)

        # Check if we can pass an array for yaw
        with context:
            vec = instrument.vectors(pos, vel, yaw=[0], **kwargs)

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


def test_local_frame_nadir_is_perpendicular_to_ellipsoid():
    """_local_frame nadir must point from satellite to geodetic sub-satellite point.

    For a satellite over 45°N: the true geodetic nadir direction is the inward
    normal to the WGS84 ellipsoid, which differs from the geocentric direction
    (pos / |pos|) by up to ~0.2° due to Earth's oblateness.  The test checks
    that the nadir vector is anti-parallel to the outward ellipsoid normal at the
    sub-satellite point, to within 0.01° (a geocentric approximation would fail
    this check at the ~0.19° level).
    """
    from pyorbital.geoloc import _local_frame

    e2 = 0.00669437999014  # WGS84 first eccentricity squared
    a = 6378.137
    phi = np.deg2rad(45.0)
    n = a / np.sqrt(1 - e2 * np.sin(phi) ** 2)
    r_surface = np.array([n * np.cos(phi), 0.0, (1 - e2) * n * np.sin(phi)])
    # Add altitude along the outward ellipsoid normal (not the geocentric radial)
    outward_normal = np.array([np.cos(phi), 0.0, np.sin(phi)])
    pos = (r_surface + outward_normal * 883.0).reshape(3, 1)
    vel = np.array([[-np.sin(phi)], [0.0], [np.cos(phi)]]) * 7.5  # km/s southward

    nadir, _, _ = _local_frame(pos, vel, nadir_convention="geodetic")

    # Outward ellipsoid normal at geodetic lat 45° is exactly (cos45, 0, sin45)
    outward_normal = np.array([np.cos(phi), 0.0, np.sin(phi)])
    dot = float(np.dot(nadir[:, 0], outward_normal))
    angle_deg = np.degrees(np.arccos(np.clip(-dot, -1, 1)))
    assert angle_deg < 0.01, f"nadir deviates {angle_deg:.4f}° from geodetic normal (should be < 0.01°)"


@pytest.mark.parametrize(("convention", "warns"), NADIR_CONFIG_CASES)
def test_arbitrary_point_geoloc(convention, warns):
    """Test geolocating an arbitrary point in the swath."""
    context = _expect_convention_warning(warns)
    with _configured_convention(convention), context:
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

        # Each convention pins its own ground solution.  The legacy and
        # geocentric answers differ by ~0.002 degrees here (about 200 m), which is
        # the effect this parametrisation exists to keep visible.
        expected = {
            # the legacy values are those released pyorbital produces: the
            # pre-0abe19f revision of this test asserted -34.69996894 / 56.69799502
            "legacy": ((-34.6999689108117, 56.6979949863614),
                       (-27.5730527370977, 55.6267408763376)),
            "geocentric": ((-34.6982153937553, 56.6970880483636),
                           (-27.5720035629012, 55.6258963697574)),
        }[convention or "legacy"]
        assert lons[0] == pytest.approx(expected[0][0])
        assert lats[0] == pytest.approx(expected[0][1])

        assert lons[2] == pytest.approx(expected[1][0])
        assert lats[2] == pytest.approx(expected[1][1])


@pytest.mark.parametrize(("convention", "warns"), NADIR_CONFIG_CASES)
def test_minimize_geoloc_error(convention, warns):
    """Test minimizing the distance to a set of gcps."""
    context = _expect_convention_warning(warns)
    # these exercise a non-zero pitch, so pin the rotation order and let the
    # parametrisation vary only the nadir convention
    with _configured_convention(convention), config.set(rotation_order="legacy"), context:
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


@pytest.mark.parametrize(("convention", "warns"), NADIR_CONFIG_CASES)
def test_minimize_time_error(convention, warns):
    """Test minimizing the distance to a set of gcps using only time offset."""
    context = _expect_convention_warning(warns)
    # these exercise a non-zero pitch, so pin the rotation order and let the
    # parametrisation vary only the nadir convention
    with _configured_convention(convention), config.set(rotation_order="legacy"), context:
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


def test_olci_scan_step_scales_time_offsets():
    """Test that OLCI scan decimation preserves elapsed observation time."""
    from pyorbital.geoloc_instrument_definitions import olci

    geom = olci(2, scan_points=[0, 3999], scan_step=100)

    assert geom._times[1, 0] == np.timedelta64(4400, "ms")


def test_slstr_nadir_scan_constant():
    """Test that SLSTR_NADIR_SCAN constant produces geometry matching legacy slstr_nadir()."""
    from pyorbital.geoloc_instrument_definitions import SLSTR_NADIR_SCAN

    legacy_geom = slstr_nadir(10)
    swath = PushbroomSwath(scanline=SLSTR_NADIR_SCAN, time_sampling=np.timedelta64(0))
    new_geom = swath.scan_geometry(scan_lines=slice(10))
    np.testing.assert_allclose(new_geom.fovs, legacy_geom.fovs)
    np.testing.assert_equal(new_geom._times, legacy_geom._times)


@pytest.mark.parametrize(("convention", "warns"), NADIR_CONFIG_CASES)
def test_bounding_box_returns_closed_polygon(convention, warns):
    """Test that bounding_box returns a closed polygon of lon/lat points."""
    context = _expect_convention_warning(warns)
    with _configured_convention(convention), context:
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


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_yaw_steering_applied_in_vectors(kwargs, warns):
    """Test that vectors() applies yaw steering when enabled."""
    context = _expect_convention_warning(warns)
    with context:
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

        vec_no_steering = instrument.vectors(pos, vel, **kwargs)
        vec_with_steering = instrument.vectors(pos, vel, yaw_steering=True, **kwargs)

        # Vectors should differ when yaw steering is applied
        assert not np.allclose(vec_no_steering, vec_with_steering)

        # Manually compute: yaw steering should add the computed yaw angle
        yaw_angle = compute_yaw_steering(pos, vel)
        vec_manual_yaw = instrument.vectors(pos, vel, yaw=yaw_angle, **kwargs)
        np.testing.assert_allclose(vec_with_steering, vec_manual_yaw)


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_compute_pixels_with_yaw_steering(kwargs, warns):
    """Test that compute_pixels passes yaw_steering through to vectors."""
    context = _expect_convention_warning(warns)
    with context:
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

        pixels_no_yaw = compute_pixels((tle1, tle2), sgeom, s_times, **kwargs)
        pixels_yaw = compute_pixels((tle1, tle2), sgeom, s_times, yaw_steering=True, **kwargs)

        # The positions should differ with yaw steering
        assert not np.allclose(pixels_no_yaw, pixels_yaw)


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_compute_pixels_with_2d_scan_times(kwargs, warns):
    """Test that compute_pixels works with multi-scan 2D time arrays (avhrr-style).

    avhrr() produces fovs of shape (2, N_SCANS, N_PX) and times of shape
    (N_SCANS, N_PX). compute_pixels must handle this without crashing and
    return pixel positions of the correct shape.
    """
    context = _expect_convention_warning(warns)
    with context:
        from pyorbital.geoloc import compute_pixels, get_lonlatalt
        from pyorbital.geoloc_instrument_definitions import avhrr

        tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
        tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
        t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)

        n_scans, n_px = 5, 10
        sgeom = avhrr(n_scans, np.arange(n_px))
        s_times = sgeom.times(t)

        assert s_times.shape == (n_scans, n_px), "times should be 2D"

        pixels = compute_pixels((tle1, tle2), sgeom, s_times, **kwargs)
        lons, lats, alts = get_lonlatalt(pixels, s_times)

        assert pixels.shape == (3, n_scans * n_px)
        assert lons.shape == (n_scans * n_px,)
        assert np.all(np.isfinite(lons))
        assert np.all(np.abs(lats) <= 90)
        # Pixels are surface points: altitude should be within terrain range (not orbital altitude)
        assert np.all(np.abs(alts) < 10000)


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_compute_pixels_2d_matches_1d_reference(kwargs, warns):
    """Test that the 2D-times path gives the same result as the 1D reference path.

    The 1D reference path (flat ScanGeometry, per-pixel SGP4) is the
    original correct behaviour.  The 2D path uses per-scan SGP4 which
    is an approximation; the difference should be sub-pixel (< 0.01 deg).
    """
    context = _expect_convention_warning(warns)
    with context:
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
        pix_ref = compute_pixels((tle1, tle2), sgeom_ref, s_times_ref, **kwargs)
        lons_ref, lats_ref, _ = get_lonlatalt(pix_ref, s_times_ref)

        # 2D path: avhrr(1, ...) with per-scan SGP4 approximation
        sgeom_2d = avhrr(1, scan_pts)
        s_times_2d = sgeom_2d.times(t)
        pix_2d = compute_pixels((tle1, tle2), sgeom_2d, s_times_2d, **kwargs)
        lons_2d, lats_2d, _ = get_lonlatalt(pix_2d, s_times_2d)

        # Per-scan SGP4 is an approximation; error < 0.01 deg for AVHRR scan rates
        np.testing.assert_allclose(lons_2d, lons_ref, atol=0.01)
        np.testing.assert_allclose(lats_2d, lats_ref, atol=0.01)


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_geolocate_multiline_scan_matches_per_pixel_orbit_propagation(kwargs, warns):
    """A compact scan must retain the orbital motion during its sweep."""
    context = _expect_convention_warning(warns)
    with context:
        from pyorbital.geoloc import ScanGeometry, geolocate

        tle = (
            "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113",
            "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875",
        )
        start = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)
        pixel_times = np.linspace(0.0, 1.48, 129)
        cross_track = np.linspace(np.deg2rad(55.0), np.deg2rad(-55.0), 129)

        compact = ScanGeometry(
            np.stack([cross_track, np.zeros_like(cross_track)])[:, np.newaxis, :],
            pixel_times[np.newaxis, :],
        )
        exact = ScanGeometry(np.stack([cross_track, np.zeros_like(cross_track)]), pixel_times)

        compact_lon, compact_lat, _ = geolocate(tle, compact, compact.times(start), **kwargs)
        exact_lon, exact_lat, _ = geolocate(tle, exact, exact.times(start), **kwargs)

        np.testing.assert_allclose(compact_lon, exact_lon, atol=9e-5)
        np.testing.assert_allclose(compact_lat, exact_lat, atol=9e-5)


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_geolocate_compact_scan_samples_only_orbit_endpoints(kwargs, warns):
    """Compact scans must not require one propagated orbit state per pixel."""
    context = _expect_convention_warning(warns)
    with context:
        from pyorbital.geoloc import ScanGeometry, geolocate

        epoch = np.datetime64("2020-01-01T00:00:00")

        class EndpointOrbit:
            def get_position(self, times, normalize=False):
                assert np.asarray(times).size <= 2
                seconds = (np.asarray(times) - epoch) / np.timedelta64(1, "s")
                position = np.vstack([np.full_like(seconds, 7000.0), 7.5 * seconds, np.zeros_like(seconds)])
                velocity = np.vstack([np.zeros_like(seconds), np.full_like(seconds, 7.5), np.zeros_like(seconds)])
                return position, velocity

        pixel_times = np.linspace(0.0, 1.48, 17)
        geometry = ScanGeometry(
            np.zeros((2, 1, pixel_times.size)),
            pixel_times[np.newaxis, :],
        )

        lon, lat, alt = geolocate(EndpointOrbit(), geometry, geometry.times(epoch), **kwargs)

        assert np.all(np.isfinite((lon, lat, alt)))


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_geolocate_compact_pass_bounds_orbit_sampling_batches(kwargs, warns):
    """A compact pass must bound temporary orbit-state batches."""
    context = _expect_convention_warning(warns)
    with context:
        from pyorbital.geoloc import ScanGeometry, geolocate

        epoch = np.datetime64("2020-01-01T00:00:00")

        class BoundedOrbit:
            def get_position(self, times, normalize=False):
                assert np.asarray(times).size <= 32
                seconds = (np.asarray(times) - epoch) / np.timedelta64(1, "s")
                position = np.vstack([np.full_like(seconds, 7000.0), 7.5 * seconds, np.zeros_like(seconds)])
                velocity = np.vstack([np.zeros_like(seconds), np.full_like(seconds, 7.5), np.zeros_like(seconds)])
                return position, velocity

        pixel_offsets = np.linspace(0.0, 1.48, 17)
        row_offsets = np.arange(100)[:, np.newaxis] * 1.5 + pixel_offsets
        geometry = ScanGeometry(np.zeros((2, 100, 17)), row_offsets)

        lon, lat, alt = geolocate(BoundedOrbit(), geometry, geometry.times(epoch), **kwargs)

        assert np.all(np.isfinite((lon, lat, alt)))


def test_get_sensor_angles_accepts_orbit_provider():
    """Sensor zenith and azimuth work with any public orbit provider."""
    from pyorbital.geoloc import get_sensor_angles
    from pyorbital.orbital import Orbital

    tle = (
        "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113",
        "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875",
    )
    time = np.datetime64("2012-12-12T04:16:01.575")
    orbit = Orbital("satellite", line1=tle[0], line2=tle[1])
    longitude = np.array([12.0])
    latitude = np.array([55.0])
    expected_azimuth, expected_elevation = orbit.get_observer_look(time, longitude, latitude, 0.0)

    zenith, azimuth = get_sensor_angles(orbit, time, longitude, latitude)

    np.testing.assert_allclose([zenith, azimuth], [90.0 - expected_elevation, expected_azimuth])


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_geolocate_accepts_per_scan_attitude(kwargs, warns):
    """Per-scan RPY arrays match independent single-scan geolocation."""
    context = _expect_convention_warning(warns)
    with context:
        from pyorbital.geoloc import geolocate
        from pyorbital.geoloc_instrument_definitions import MultiLineWhiskbroomScan

        tle = (
            "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113",
            "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875",
        )
        start = np.datetime64("2012-12-12T04:16:01.575")
        scan = MultiLineWhiskbroomScan(20, 55.0, 1.5, 0.001, 1, 1.0 / 830.0)
        attitudes = np.array([[0.001, 0.0, 0.0], [-0.001, 0.0, 0.0]])
        combined_geometry = scan.scan_geometry(2)

        combined = geolocate(tle, combined_geometry, combined_geometry.times(start), rpy=attitudes, **kwargs)
        independent = [
            geolocate(tle, geometry, geometry.times(start), rpy=attitude, **kwargs)
            for geometry, attitude in zip(
                [scan.scan_geometry(1, scan_offsets=[0.0]), scan.scan_geometry(1, scan_offsets=[1.5])],
                attitudes,
                strict=True,
            )
        ]

        np.testing.assert_allclose(combined[:2], [np.concatenate([part[i] for part in independent]) for i in range(2)])


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
    expected = ((lines_per_scan - 1) / 2.0 - np.arange(lines_per_scan)) * step
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

@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_compute_pixel_works(kwargs, warns):
    """Check that compute_pixels works for MultiLineWhiskbroomScan without OOM or crash."""
    context = _expect_convention_warning(warns)
    # MERSI carries non-zero along-track detector angles, so the rotation order
    # is in play here too; pin it and vary only the nadir convention
    with config.set(rotation_order="legacy"), context:
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
        pixels = compute_pixels((tle1, tle2), sgeom, s_times, **kwargs)
        lons, lats, alts = get_lonlatalt(pixels, s_times)

        assert pixels.shape == (3, n_scans * 40 * 8192)
        assert np.all(np.isfinite(lons))
        assert np.all(np.abs(lats) <= 90)


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_geolocate_matches_compute_pixels_get_lonlatalt(kwargs, warns):
    """Check that short-arc geolocation agrees with per-pixel propagation within 10 m."""
    context = _expect_convention_warning(warns)
    # these exercise a non-zero pitch, so pin the rotation order and let the
    # parametrisation vary only the nadir convention
    with config.set(rotation_order="legacy"), context:
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

        # Reference: full orbit propagation at every pixel acquisition time
        flat_geometry = ScanGeometry(sgeom.fovs.reshape(2, -1), s_times.reshape(-1) - s_times.reshape(-1)[0])
        pixels = compute_pixels((tle1, tle2), flat_geometry, s_times.reshape(-1), **kwargs)
        lon_ref, lat_ref, alt_ref = get_lonlatalt(pixels, s_times.reshape(-1))

        # New fused path
        lon_fused, lat_fused, alt_fused = geolocate((tle1, tle2), sgeom, s_times, **kwargs)

        assert lon_fused.shape == lon_ref.shape
        np.testing.assert_allclose(lon_fused, lon_ref, atol=9e-5,
                                   err_msg="geolocate lon disagrees with reference")
        np.testing.assert_allclose(lat_fused, lat_ref, atol=1e-6,
                                   err_msg="geolocate lat disagrees with reference")
        np.testing.assert_allclose(alt_fused, alt_ref, atol=1.0,
                                   err_msg="geolocate alt disagrees with reference")


@pytest.mark.parametrize(("kwargs", "warns"), NADIR_CONVENTION_CASES)
def test_geolocate_with_yaw_steering_matches_reference(kwargs, warns):
    """Check yaw-steered short-arc geolocation against per-pixel propagation."""
    context = _expect_convention_warning(warns)
    # these exercise a non-zero pitch, so pin the rotation order and let the
    # parametrisation vary only the nadir convention
    with config.set(rotation_order="legacy"), context:
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

        # Reference: full orbit propagation at every pixel acquisition time
        flat_geometry = ScanGeometry(sgeom.fovs.reshape(2, -1), s_times.reshape(-1) - s_times.reshape(-1)[0])
        pixels = compute_pixels((tle1, tle2), flat_geometry, s_times.reshape(-1), yaw_steering=True, **kwargs)
        lon_ref, lat_ref, alt_ref = get_lonlatalt(pixels, s_times.reshape(-1))

        # Fused numba path with yaw steering
        lon_f, lat_f, alt_f = geolocate((tle1, tle2), sgeom, s_times, yaw_steering=True, **kwargs)

        assert lon_f.shape == lon_ref.shape
        np.testing.assert_allclose(lon_f, lon_ref, atol=9e-5,
                                   err_msg="yaw-steered lon disagrees with reference")
        np.testing.assert_allclose(lat_f, lat_ref, atol=1e-6,
                                   err_msg="yaw-steered lat disagrees with reference")
        np.testing.assert_allclose(alt_f, alt_ref, atol=1.0,
                                   err_msg="yaw-steered alt disagrees with reference")


# ---------------------------------------------------------------------------
# FocalPlaneWhiskbroomScan tests
# ---------------------------------------------------------------------------

def _fp_scan(lines_per_scan=4, altitude_km=830.0, det_position=(0.0, 0.0), boresight=None):
    """Helper: FocalPlaneWhiskbroomScan with focal_length=altitude_km, det_space=1.0."""
    return FocalPlaneWhiskbroomScan(
        pixels_per_scan=20,
        scan_angle=55.0,
        scan_rate=1.5,
        pixel_dwell_time=0.001,
        lines_per_scan=lines_per_scan,
        focal_length=altitude_km,
        det_space=1.0,
        det_position=det_position,
        boresight=boresight,
    )


def test_focal_plane_identity_matches_multiline_whiskbroom():
    """FocalPlaneWhiskbroomScan with identity boresight and zero offset matches MultiLineWhiskbroomScan."""
    altitude_km = 830.0
    lines_per_scan = 4
    fp = _fp_scan(lines_per_scan=lines_per_scan, altitude_km=altitude_km)
    ref = MultiLineWhiskbroomScan(
        pixels_per_scan=20,
        scan_angle=55.0,
        scan_rate=1.5,
        pixel_dwell_time=0.001,
        lines_per_scan=lines_per_scan,
        along_track_step=1.0 / altitude_km,
    )
    # FocalPlaneWhiskbroomScan uses exact arctan2; MultiLineWhiskbroomScan uses the
    # small-angle approximation tan(θ)≈θ.  Tolerance covers that ~nrad difference.
    np.testing.assert_allclose(fp.along_track_angles(), ref.along_track_angles(), atol=1e-7)
    np.testing.assert_allclose(fp.cross_track_angles(), ref.cross_track_angles(), atol=1e-10)


def test_focal_plane_det_position_y_shifts_along_track():
    """det_position y-offset shifts each row's along-track angle by the correct arctan amount."""
    altitude_km = 830.0
    y_offset_mm = 5.0
    lines_per_scan = 4
    fp_base = _fp_scan(lines_per_scan=lines_per_scan, altitude_km=altitude_km, det_position=(0.0, 0.0))
    fp_off = _fp_scan(lines_per_scan=lines_per_scan, altitude_km=altitude_km, det_position=(0.0, y_offset_mm))

    # Compute expected angles from first principles
    rows = np.arange(lines_per_scan, dtype=float)
    y_base = ((lines_per_scan - 1) / 2.0 - rows) * 1.0          # det_space=1.0
    y_off = y_base + y_offset_mm
    expected_base = np.arctan2(y_base, altitude_km)
    expected_off = np.arctan2(y_off, altitude_km)

    np.testing.assert_allclose(fp_base.along_track_angles(), expected_base, atol=1e-12)
    np.testing.assert_allclose(fp_off.along_track_angles(), expected_off, atol=1e-12)


def test_focal_plane_accepts_effective_detector_along_track_offsets():
    """Effective detector offsets augment nominal focal-plane angles."""
    from dataclasses import replace

    base = _fp_scan(lines_per_scan=4)
    offsets = np.array([[0.0, -0.003], [0.0, -0.001], [0.0, 0.001], [0.0, 0.003]])
    corrected = replace(base, detector_offsets_rad=offsets)

    np.testing.assert_allclose(corrected.along_track_angles(), base.along_track_angles() + offsets[:, 1])
    np.testing.assert_allclose(corrected.scan_geometry(1).fovs[1], base.scan_geometry(1).fovs[1] + offsets[:, 1:])


def test_focal_plane_accepts_effective_detector_cross_track_offsets():
    """Cross-track detector offsets are retained in the generated geometry."""
    from dataclasses import replace

    base = _fp_scan(lines_per_scan=4)
    offsets = np.array([[-0.003, 0.0], [-0.001, 0.0], [0.001, 0.0], [0.003, 0.0]])
    corrected = replace(base, detector_offsets_rad=offsets)
    expected = base.cross_track_angles()[np.newaxis, :] + offsets[:, :1]

    np.testing.assert_allclose(corrected.scan_geometry(1).fovs[0], expected)


def test_focal_plane_accepts_polynomial_scan_angle_corrections():
    """Odd and even scan-law corrections use normalized scan position."""
    from dataclasses import replace

    base = _fp_scan()
    coefficients = (0.0, 0.002, 0.003)
    corrected = replace(base, scan_angle_correction_coefficients_rad=coefficients)
    scan_position = np.linspace(1.0, -1.0, base.pixels_per_scan)
    expected = base.cross_track_angles() + np.polynomial.polynomial.polyval(scan_position, coefficients)

    np.testing.assert_allclose(corrected.cross_track_angles(), expected)


def test_focal_plane_accepts_detector_scan_position_slopes():
    """Detector LOS slopes vary linearly across normalized scan position."""
    from dataclasses import replace

    base = _fp_scan(lines_per_scan=4)
    slopes = np.array([[0.001, -0.002], [0.002, -0.001], [-0.001, 0.002], [-0.002, 0.001]])
    corrected = replace(base, detector_scan_slopes_rad=slopes)
    scan_position = np.linspace(1.0, -1.0, base.pixels_per_scan)
    expected = base.scan_geometry(1).fovs + np.stack(
        [slopes[:, :1] * scan_position, slopes[:, 1:] * scan_position],
    )

    np.testing.assert_allclose(corrected.scan_geometry(1).fovs, expected)


def test_focal_plane_scan_geometry_accepts_explicit_scan_offsets():
    """Telemetry-derived scan offsets replace nominal scan-rate spacing."""
    scan = _fp_scan(lines_per_scan=2)

    geometry = scan.scan_geometry(2, scan_offsets=np.array([0.25, 2.25]))

    np.testing.assert_allclose(
        geometry._times.astype("timedelta64[ns]").astype(float)[:, 0] / 1e9,
        [0.25, 0.25, 2.25, 2.25],
    )


def test_focal_plane_scan_geometry_accepts_per_scan_los_offsets():
    """Per-scan LOS offsets can represent mirror-side geometry."""
    scan = _fp_scan(lines_per_scan=2)
    los_offsets = np.array([[0.001, -0.002], [-0.003, 0.004]])
    base = scan.scan_geometry(2).fovs.reshape(2, 2, 2, scan.pixels_per_scan)

    corrected = scan.scan_geometry(2, scan_los_offsets_rad=los_offsets).fovs.reshape(base.shape)

    np.testing.assert_allclose(corrected - base, np.broadcast_to(los_offsets.T[:, :, None, None], base.shape))


def test_focal_plane_det_position_x_shifts_cross_track_uniformly():
    """det_position x-offset shifts all cross-track angles by the same amount."""
    altitude_km = 830.0
    x_offset_mm = 3.0
    fp_base = _fp_scan(altitude_km=altitude_km, det_position=(0.0, 0.0))
    fp_off = _fp_scan(altitude_km=altitude_km, det_position=(x_offset_mm, 0.0))

    delta = fp_off.cross_track_angles() - fp_base.cross_track_angles()
    np.testing.assert_allclose(delta, delta[0], atol=1e-12,
                               err_msg="cross-track shift is not uniform across pixels")
    expected = np.arctan2(x_offset_mm, altitude_km)
    np.testing.assert_allclose(delta[0], expected, atol=1e-10,
                               err_msg="cross-track shift magnitude is incorrect")


def test_focal_plane_pure_pitch_boresight_shifts_along_track():
    """A pure pitch rotation in the boresight shifts all along-track angles by that angle."""
    pitch_rad = np.deg2rad(0.05)
    R_pitch = np.array([
        [1.0, 0.0,            0.0           ],
        [0.0, np.cos(pitch_rad), -np.sin(pitch_rad)],
        [0.0, np.sin(pitch_rad),  np.cos(pitch_rad)],
    ])
    fp_no = _fp_scan(boresight=None)
    fp_bs = _fp_scan(boresight=R_pitch)

    delta = fp_bs.along_track_angles() - fp_no.along_track_angles()
    # Rx(+θ) rotates the nadir direction toward −y, so the apparent along-track
    # angle of the boresight shifts by −θ (backward).
    np.testing.assert_allclose(delta, -pitch_rad, atol=1e-6,
                               err_msg="pitch boresight does not shift along-track angles correctly")


def test_focal_plane_scan_geometry_shape_and_times():
    """scan_geometry produces correctly shaped fovs and times arrays."""
    n_scans, L, N = 3, 4, 20
    fp = _fp_scan(lines_per_scan=L)
    geom = fp.scan_geometry(n_scans=n_scans)
    assert geom.fovs.shape == (2, n_scans * L, N)
    assert geom._times.shape == (n_scans * L, N)
    assert geom.lines_per_scan == L


def test_focal_plane_structural_invariants_preserved():
    """Cross-track angles are identical for all rows; along-track is constant across pixels."""
    fp = _fp_scan(lines_per_scan=4, det_position=(1.5, 2.0))
    geom = fp.scan_geometry(n_scans=2)
    # Cross-track: all rows in a scan must be the same
    for row in range(1, 4):
        np.testing.assert_array_equal(geom.fovs[0, 0, :], geom.fovs[0, row, :],
                                      err_msg=f"cross-track row {row} differs from row 0")
    # Along-track: constant across pixels within each row
    along = geom.fovs[1]         # (n_scans*L, N)
    assert np.all(along == along[:, 0:1]), "along-track angles vary across pixels"


@pytest.mark.parametrize("rotation_order", ["legacy", "pitch_first"])
def test_along_track_spacing_depends_on_the_rotation_order(rotation_order):
    """Along-track pixel spacing must not depend on the cross-track scan angle.

    Before the rotation-order fix, the along-track component of the pointing
    vector was compressed by cos(scan_angle), giving ~57% of the correct value
    at 55°.  After the fix, along_track · r == -sin(alpha) regardless of theta.

    Geometry: equatorial satellite above the x-axis with northward velocity.
    This gives a clean analytical frame:
        nadir       = [-1,  0,  0]
        along_track = [ 0,  1,  0]  (northward)
        cross_track = [ 0,  0, -1]
    so the expected along-track component is exactly -sin(alpha).

    Uses the 3-D broadcast path (fovs.shape = (2, n_rows, n_pixels)) so that
    along-track angles vary per row while cross-track angles vary per pixel.
    """
    with config.set(nadir_convention="geocentric"):
        alpha = 1.0 / 830.0        # typical MERSI along-track step (rad)
        theta = np.deg2rad(55.13)  # full-swath edge scan angle

        # 2 rows x 2 pixels:
        #   row 0: no along-track offset;  row 1: alpha offset
        #   col 0: nadir (theta=0);        col 1: swath edge (theta)
        fovs = np.array([
            [[0.0, theta], [0.0, theta]],   # fovs[0]: cross-track, same per row
            [[0.0, 0.0],   [alpha, alpha]], # fovs[1]: along-track, same per col
        ])  # shape (2, 2, 2)
        geom = ScanGeometry(fovs, np.zeros((2, 2)), lines_per_scan=2)

        # Equatorial satellite with northward velocity; use same position for both rows
        # (rows are separated by microseconds — no measurable orbital movement).
        R_orbit = 7208.0  # km (altitude ~830 km)
        pos = np.array([[R_orbit, R_orbit], [0.0, 0.0], [0.0, 0.0]])
        vel = np.array([[0.0, 0.0], [7.34, 7.34], [0.0, 0.0]])  # km/s, northward

        vectors = geom.vectors(pos, vel,
                               rotation_order=rotation_order)  # (3, 4): [r0-nadir, r0-edge, r1-nadir, r1-edge]

        # For this geometry, along_track = [0, 1, 0] analytically.
        along_track = np.array([0.0, 1.0, 0.0])

        at_r0_nadir = along_track @ vectors[:, 0]  # row 0, nadir:  expect ~0
        at_r1_nadir = along_track @ vectors[:, 2]  # row 1, nadir:  expect ~-sin(alpha)
        at_r1_edge  = along_track @ vectors[:, 3]  # row 1, edge:   expect ~-sin(alpha)

        np.testing.assert_allclose(at_r0_nadir, 0.0,           atol=1e-6,
                                   err_msg="row-0 nadir along-track component should be ~0")
        np.testing.assert_allclose(at_r1_nadir, -np.sin(alpha), rtol=1e-5,
                                   err_msg="row-1 nadir along-track component should be -sin(alpha)")
        # This is exactly where the two orders part company: applying the pitch
        # first leaves the along-track component untouched by the sweep, while
        # applying the roll first compresses it by cos(theta) -- 57% of the value
        # at 55 degrees, which is what the later order was introduced to avoid.
        expected_edge = (-np.sin(alpha) if rotation_order == "pitch_first"
                         else -np.cos(theta) * np.sin(alpha))
        np.testing.assert_allclose(at_r1_edge, expected_edge, rtol=1e-5,
                                   err_msg=f"swath-edge along-track component, {rotation_order}")


def test_local_frame_geocentric_nadir_points_at_earth_centre():
    """The geocentric convention must give a nadir exactly anti-parallel to the position.

    Validation against reference geolocation shows FY-3 references its attitude to
    the geocentric (orbital) local vertical rather than the ellipsoid normal, so
    this convention has to be available and has to be exact: any admixture of the
    geodetic normal reintroduces an oblateness-driven error of order a kilometre
    on the ground.
    """
    from pyorbital.geoloc import _local_frame

    e2 = 0.00669437999014
    a = 6378.137
    phi = np.deg2rad(45.0)
    n = a / np.sqrt(1 - e2 * np.sin(phi) ** 2)
    outward_normal = np.array([np.cos(phi), 0.0, np.sin(phi)])
    pos = (np.array([n * np.cos(phi), 0.0, (1 - e2) * n * np.sin(phi)]) + outward_normal * 883.0).reshape(3, 1)
    # a circular-orbit velocity, i.e. perpendicular to the position: the frame
    # re-orthogonalises the nadir against the velocity, so a velocity that is not
    # perpendicular to pos would legitimately tilt the nadir away from the centre
    # A realistic inclined-orbit velocity: perpendicular to the position but NOT
    # coplanar with the meridian.  With a coplanar velocity the Gram-Schmidt step
    # projects every nadir convention onto -pos/|pos| and they become
    # indistinguishable, so a polar-plane test geometry proves nothing.
    orbit_normal = np.array([0.0, -np.sin(np.deg2rad(98.7)), np.cos(np.deg2rad(98.7))])
    vel = np.cross(orbit_normal, pos[:, 0]).reshape(3, 1)
    vel = vel / np.sqrt((vel**2).sum()) * 7.5

    nadir, _, _ = _local_frame(pos, vel, nadir_convention="geocentric")

    expected = -pos[:, 0] / np.sqrt((pos[:, 0] ** 2).sum())
    np.testing.assert_allclose(nadir[:, 0], expected, rtol=1e-12, atol=1e-12)


def test_local_frame_defaults_to_the_released_nadir_convention():
    """The default must stay the historical convention so existing products reproduce.

    Released pyorbital builds the nadir as the normalised position of the
    ellipsoid subpoint of the antipode.  Archives such as the reprocessed AVHRR
    record were produced with it, and their fitted attitude corrections absorbed
    whatever it implied, so changing the default silently would move every one of
    those products.  The corrected conventions are opt-in until a major release.
    """
    from pyorbital.geoloc import _local_frame

    e2 = 0.00669437999014
    a = 6378.137
    phi = np.deg2rad(45.0)
    n = a / np.sqrt(1 - e2 * np.sin(phi) ** 2)
    outward_normal = np.array([np.cos(phi), 0.0, np.sin(phi)])
    pos = (np.array([n * np.cos(phi), 0.0, (1 - e2) * n * np.sin(phi)]) + outward_normal * 883.0).reshape(3, 1)
    # A realistic inclined-orbit velocity: perpendicular to the position but NOT
    # coplanar with the meridian.  With a coplanar velocity the Gram-Schmidt step
    # projects every nadir convention onto -pos/|pos| and they become
    # indistinguishable, so a polar-plane test geometry proves nothing.
    orbit_normal = np.array([0.0, -np.sin(np.deg2rad(98.7)), np.cos(np.deg2rad(98.7))])
    vel = np.cross(orbit_normal, pos[:, 0]).reshape(3, 1)
    vel = vel / np.sqrt((vel**2).sum()) * 7.5

    with _expect_convention_warning(warns=True):
        default_nadir, _, _ = _local_frame(pos, vel)
    legacy_nadir, _, _ = _local_frame(pos, vel, nadir_convention="legacy")
    geocentric_nadir, _, _ = _local_frame(pos, vel, nadir_convention="geocentric")

    np.testing.assert_allclose(default_nadir, legacy_nadir, rtol=1e-12, atol=1e-12)
    # and it is genuinely a different convention, not an alias of the corrected one
    assert not np.allclose(default_nadir, geocentric_nadir, atol=1e-9)


def test_local_frame_explicit_legacy_choice_does_not_warn():
    """Asking for the legacy convention on purpose is an informed choice, not a mistake.

    The DeprecationWarning exists to catch callers who are silently relying on the
    default.  Someone reprocessing an existing archive has to pin the old
    convention deliberately, and warning at them on every call would be noise --
    especially as the suite runs with ``filterwarnings = error``.
    """
    from pyorbital.geoloc import _local_frame

    pos = np.array([[4507.0], [0.0], [4497.0]])
    vel = np.array([[-4497.0], [0.0], [4507.0]])
    vel = vel / np.sqrt((vel**2).sum()) * 7.5

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _local_frame(pos, vel, nadir_convention="legacy")


# Each nadir convention, with whether relying on it raises the deprecation.
# Parametrising over these keeps every geometry test honest about which
@pytest.mark.parametrize(("kwargs", "warns", "convention"), NADIR_CONVENTIONS)
def test_scan_geometry_vectors_honours_the_nadir_convention(kwargs, warns, convention):
    """``ScanGeometry.vectors`` must forward the convention, not just ``_local_frame``.

    Callers reach the geometry through the public entry points, so the choice has
    to be selectable there; otherwise the only way to opt into the corrected
    computation is to monkeypatch an underscore-private function.
    """
    from pyorbital.geoloc import _local_frame

    scan_angles = np.deg2rad(np.array([-55.0, 0.0, 55.0]))
    geometry = ScanGeometry(
        np.vstack((scan_angles, np.zeros_like(scan_angles))),
        np.array([0.0, 0.1, 0.2]),
    )
    pos, vel = _inclined_orbit_state()
    # vectors() works per pixel, so the orbit state is repeated across the scan
    pos = np.repeat(pos, scan_angles.size, axis=1)
    vel = np.repeat(vel, scan_angles.size, axis=1)

    context = _expect_convention_warning(warns)
    with context:
        vectors = geometry.vectors(pos, vel, **kwargs)

    # the centre pixel is at zero scan angle, so it must lie along the nadir of
    # whichever convention was selected
    expected_nadir, _, _ = _local_frame(pos, vel, nadir_convention=convention)
    np.testing.assert_allclose(vectors[:, 1], expected_nadir[:, 1], rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize(("kwargs", "warns", "convention"), NADIR_CONVENTIONS)
def test_compute_pixels_honours_the_nadir_convention(kwargs, warns, convention):
    """``compute_pixels`` must expose the convention to its callers.

    pygac and other downstream users reach the geometry through this function,
    so it is the level at which an archive reprocessing pins its choice.
    """
    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
    scan_angles = np.deg2rad(np.array([-55.0, 0.0, 55.0]))
    geometry = ScanGeometry(
        np.vstack((scan_angles, np.zeros_like(scan_angles))),
        np.array([0.0, 0.1, 0.2]),
    )
    times = geometry.times(dt.datetime(2012, 12, 12, 4, 16, 1, 575000))

    context = _expect_convention_warning(warns)
    with context:
        pixels = compute_pixels((tle1, tle2), geometry, times, **kwargs)

    reference = compute_pixels((tle1, tle2), geometry, times, nadir_convention=convention)
    np.testing.assert_allclose(pixels, reference, rtol=1e-12, atol=1e-12)
    if convention == "geocentric":
        legacy = compute_pixels((tle1, tle2), geometry, times, nadir_convention="legacy")
        assert not np.allclose(pixels, legacy, atol=1e-6)


@pytest.mark.parametrize(("kwargs", "warns", "convention"), NADIR_CONVENTIONS)
def test_geolocate_honours_the_nadir_convention(kwargs, warns, convention):
    """``geolocate`` is the top-level entry point and must carry the choice through.

    It dispatches to several internal paths (fused, per-pixel-orbit, per-scan
    attitude); the convention has to survive whichever one is taken.
    """
    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
    scan_angles = np.deg2rad(np.array([-55.0, 0.0, 55.0]))
    geometry = ScanGeometry(
        np.vstack((scan_angles, np.zeros_like(scan_angles))),
        np.array([0.0, 0.1, 0.2]),
    )
    times = geometry.times(dt.datetime(2012, 12, 12, 4, 16, 1, 575000))

    context = _expect_convention_warning(warns)
    with context:
        lon, lat, alt = geolocate(Orbital("mysatellite", line1=tle1, line2=tle2),
                                  geometry, times, **kwargs)

    ref_lon, ref_lat, _ = geolocate(Orbital("mysatellite", line1=tle1, line2=tle2),
                                    geometry, times, nadir_convention=convention)
    np.testing.assert_allclose(lon, ref_lon, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(lat, ref_lat, rtol=1e-12, atol=1e-12)


def test_nadir_convention_can_be_set_through_the_config():
    """The convention must be selectable without touching every call site.

    Downstream libraries (pygac and the readers built on it) call into the
    geolocation without forwarding a convention argument, so an archive
    reprocessing has to be able to choose one process-wide.  Setting it through
    the config is a deliberate choice, so it must not warn.
    """
    from pyorbital.config import config
    from pyorbital.geoloc import _local_frame

    pos, vel = _inclined_orbit_state()

    with config.set(nadir_convention="geocentric"), warnings.catch_warnings():
        warnings.simplefilter("error")
        configured, _, _ = _local_frame(pos, vel)

    expected, _, _ = _local_frame(pos, vel, nadir_convention="geocentric")
    np.testing.assert_allclose(configured, expected, rtol=1e-12, atol=1e-12)


def test_explicit_nadir_convention_overrides_the_config():
    """An argument at the call site beats the process-wide setting.

    Comparing conventions against reference geolocation means running both in one
    process, so the explicit form has to win over whatever the config says.
    """
    from pyorbital.config import config
    from pyorbital.geoloc import _local_frame

    pos, vel = _inclined_orbit_state()

    with config.set(nadir_convention="geocentric"):
        overridden, _, _ = _local_frame(pos, vel, nadir_convention="legacy")

    expected, _, _ = _local_frame(pos, vel, nadir_convention="legacy")
    np.testing.assert_allclose(overridden, expected, rtol=1e-12, atol=1e-12)


ROTATION_ORDER_CASES = [
    pytest.param({}, True, id="default-warns"),
    pytest.param({"rotation_order": "legacy"}, False, id="explicit-legacy"),
    pytest.param({"rotation_order": "pitch_first"}, False, id="pitch-first"),
]


def _expect_rotation_warning(warns):
    """Return a context asserting the legacy rotation-order deprecation fires, or not."""
    if warns:
        return pytest.warns(DeprecationWarning, match="legacy rotation order")
    return contextlib.nullcontext()


@pytest.mark.parametrize(("kwargs", "warns"), ROTATION_ORDER_CASES)
def test_scan_geometry_vectors_honours_the_rotation_order(kwargs, warns):
    """Roll and pitch do not commute, so the order has to be selectable.

    The cross-track fov carries the whole scan angle while pitch is a small
    attitude bias, so the two orders diverge with scan angle: nothing at nadir,
    growing to ~1.8 km at the AVHRR swath edge for a 0.16 degree pitch.  The
    released order applies roll first; commit 0abe19f swapped it.
    """
    scan_angles = np.deg2rad(np.array([-55.0, 0.0, 55.0]))
    geometry = ScanGeometry(
        np.vstack((scan_angles, np.zeros_like(scan_angles))),
        np.array([0.0, 0.1, 0.2]),
    )
    pos, vel = _inclined_orbit_state()
    pos = np.repeat(pos, scan_angles.size, axis=1)
    vel = np.repeat(vel, scan_angles.size, axis=1)
    pitch = np.deg2rad(0.16)

    with _expect_rotation_warning(warns):
        vectors = geometry.vectors(pos, vel, 0.0, pitch, nadir_convention="geocentric", **kwargs)

    reference = geometry.vectors(pos, vel, 0.0, pitch, nadir_convention="geocentric",
                                 rotation_order=kwargs.get("rotation_order", "legacy"))
    np.testing.assert_allclose(vectors, reference, rtol=1e-12, atol=1e-12)

    # the two orders must genuinely differ, and only away from nadir
    legacy = geometry.vectors(pos, vel, 0.0, pitch, nadir_convention="geocentric",
                              rotation_order="legacy")
    pitch_first = geometry.vectors(pos, vel, 0.0, pitch, nadir_convention="geocentric",
                                   rotation_order="pitch_first")
    np.testing.assert_allclose(legacy[:, 1], pitch_first[:, 1], rtol=1e-9, atol=1e-9)
    assert not np.allclose(legacy[:, 0], pitch_first[:, 0], atol=1e-9)


@pytest.mark.parametrize("rotation_order", ["legacy", "pitch_first"])
def test_fused_path_honours_the_rotation_order(rotation_order):
    """The numba kernel must apply the same rotation order as the array path.

    ``geolocate`` dispatches to a fused kernel for 3-D fovs whose times are
    constant along a row (pushbroom-style geometries).  The kernel reimplements
    the rotations in closed form, so it has to be kept in step with
    ``ScanGeometry.vectors`` -- otherwise selecting an order would silently
    change the answer only on some instruments.
    """
    pytest.importorskip("numba")
    from pyorbital.geoloc import _HAS_NUMBA, compute_pixels, geolocate, get_lonlatalt

    if not _HAS_NUMBA:
        pytest.skip("numba not available")

    tle1 = "1 33591U 09005A   12345.45213434  .00000391  00000-0  24004-3 0  6113"
    tle2 = "2 33591 098.8821 283.2036 0013384 242.4835 117.4960 14.11432063197875"
    t = dt.datetime(2012, 12, 12, 4, 16, 1, 575000)

    rows, pixels_per_row = 4, 9
    cross = np.deg2rad(np.linspace(-50.0, 50.0, pixels_per_row))
    # non-zero along-track angles per row: this is what makes the order matter
    along = np.deg2rad(np.linspace(-0.6, 0.6, rows))
    fovs = np.empty((2, rows, pixels_per_row))
    fovs[0] = cross[np.newaxis, :]
    fovs[1] = along[:, np.newaxis]
    # times constant along each row, so the fused path is taken
    times = np.repeat(np.arange(rows) * 0.25, pixels_per_row).reshape(rows, pixels_per_row)
    geometry = ScanGeometry(fovs, times, lines_per_scan=1)
    scan_times = geometry.times(t)

    fused_lon, fused_lat, _ = geolocate((tle1, tle2), geometry, scan_times,
                                        nadir_convention="geocentric",
                                        rotation_order=rotation_order)

    flat = ScanGeometry(fovs.reshape(2, -1), scan_times.reshape(-1) - scan_times.reshape(-1)[0])
    reference = compute_pixels((tle1, tle2), flat, scan_times.reshape(-1),
                               nadir_convention="geocentric", rotation_order=rotation_order)
    ref_lon, ref_lat, _ = get_lonlatalt(reference, scan_times.reshape(-1))

    np.testing.assert_allclose(fused_lon, ref_lon, atol=1e-6)
    np.testing.assert_allclose(fused_lat, ref_lat, atol=1e-6)


def test_legacy_nadir_is_not_orthogonalised():
    """The legacy convention must reproduce released pyorbital exactly.

    Released pyorbital builds the triad without re-orthogonalising the nadir
    against the velocity.  The Gram-Schmidt step added later changes the ground
    solution by up to ~2.7 km once a pitch bias is involved -- more than either
    the nadir convention or the rotation order -- so "legacy" only means
    "as released" if the orthogonalisation is skipped too.
    """
    from pyorbital.geoloc import _local_frame, subpoint, vnorm

    pos, vel = _inclined_orbit_state()

    nadir, along, _ = _local_frame(pos, vel, nadir_convention="legacy")

    expected = subpoint(-pos)
    expected = expected / vnorm(expected)
    np.testing.assert_allclose(nadir, expected, rtol=1e-12, atol=1e-12)
    # released behaviour leaves the triad slightly non-orthogonal; that is the point
    assert abs(float(np.sum(nadir * along))) > 1e-6


def test_corrected_conventions_keep_an_orthonormal_frame():
    """The corrected conventions must stay orthonormal.

    The broadcast rotation relies on nadir being perpendicular to along_track,
    and rotating about a non-orthonormal triad is not well defined, so only the
    legacy reproduction path is allowed to skip it.
    """
    from pyorbital.geoloc import _local_frame

    pos, vel = _inclined_orbit_state()

    for convention in ("geocentric", "geodetic"):
        nadir, along, cross = _local_frame(pos, vel, nadir_convention=convention)
        assert abs(float(np.sum(nadir * along))) < 1e-12, convention
        assert abs(float(np.sum(nadir * cross))) < 1e-12, convention


def test_legacy_frame_does_not_use_the_orthonormal_broadcast_path():
    """The broadcast rotation assumes nadir ⊥ along_track, which legacy breaks.

    ``_vectors_broadcast`` exploits that orthogonality to simplify the
    Rodrigues rotation, so feeding it the deliberately non-orthogonal legacy
    frame would silently produce a different answer from the general path.
    """
    scan_angles = np.deg2rad(np.linspace(-50.0, 50.0, 7))
    rows = 3
    fovs = np.empty((2, rows, scan_angles.size))
    fovs[0] = scan_angles[np.newaxis, :]
    fovs[1] = np.deg2rad(np.linspace(-0.5, 0.5, rows))[:, np.newaxis]
    times = np.repeat(np.arange(rows) * 0.25, scan_angles.size).reshape(rows, scan_angles.size)
    geometry = ScanGeometry(fovs, times, lines_per_scan=1)

    pos, vel = _inclined_orbit_state()
    pos = np.repeat(pos, rows, axis=1)
    vel = np.repeat(vel, rows, axis=1)

    broadcast = geometry.vectors(pos, vel, nadir_convention="legacy", rotation_order="legacy")

    flat = ScanGeometry(fovs.reshape(2, -1), times.reshape(-1))
    per_pixel_pos = np.repeat(pos, scan_angles.size, axis=1)
    per_pixel_vel = np.repeat(vel, scan_angles.size, axis=1)
    reference = flat.vectors(per_pixel_pos, per_pixel_vel,
                             nadir_convention="legacy", rotation_order="legacy")

    np.testing.assert_allclose(broadcast, reference, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("rotation_order", ["legacy", "pitch_first"])
def test_broadcast_path_honours_the_rotation_order(rotation_order):
    """The broadcast rotation must apply the requested order like the other paths.

    It reimplements the composition in closed form for speed, so without this it
    silently applies one order regardless of what the caller asked for -- and it
    is the path taken by any 3-D scan geometry with per-scan orbit state.
    """
    scan_angles = np.deg2rad(np.linspace(-50.0, 50.0, 7))
    rows = 3
    fovs = np.empty((2, rows, scan_angles.size))
    fovs[0] = scan_angles[np.newaxis, :]
    fovs[1] = np.deg2rad(np.linspace(-0.5, 0.5, rows))[:, np.newaxis]
    times = np.repeat(np.arange(rows) * 0.25, scan_angles.size).reshape(rows, scan_angles.size)
    geometry = ScanGeometry(fovs, times, lines_per_scan=1)

    pos, vel = _inclined_orbit_state()
    pos = np.repeat(pos, rows, axis=1)
    vel = np.repeat(vel, rows, axis=1)

    broadcast = geometry.vectors(pos, vel, nadir_convention="geocentric",
                                 rotation_order=rotation_order)

    flat = ScanGeometry(fovs.reshape(2, -1), times.reshape(-1))
    reference = flat.vectors(np.repeat(pos, scan_angles.size, axis=1),
                             np.repeat(vel, scan_angles.size, axis=1),
                             nadir_convention="geocentric", rotation_order=rotation_order)

    np.testing.assert_allclose(broadcast, reference, rtol=1e-10, atol=1e-10)
