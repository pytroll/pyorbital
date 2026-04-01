"""Module to compute geolocalization of a satellite scene."""


from __future__ import print_function

import numpy as np
from pyproj import Transformer

# DIRTY STUFF. Needed the get_lonlatalt function to work on pos directly if
# we want to print out lonlats in the end.
from pyorbital import astronomy
from pyorbital.orbital import Orbital

A = 6378.137  # WGS84 Equatorial radius (km)
B = 6356.75231414  # km, GRS80

OMEGA_EARTH = 7.2921159e-5  # Earth's rotation rate (rad/s)


def geodetic_lat(point, a=A, b=B):
    """Get the Geodetic latitude of a point."""
    x, y, z = point
    r = np.sqrt(x * x + y * y)
    geoc_lat = np.arctan2(z, r)

    geod_lat = geoc_lat
    e2 = (a * a - b * b) / (a * a)
    while True:
        phi = geod_lat
        C = 1 / np.sqrt(1 - e2 * np.sin(phi)**2)
        geod_lat = np.arctan2(z + a * C * e2 * np.sin(phi), r)
        if np.allclose(geod_lat, phi):
            return geod_lat


def subpoint(query_point, a=A, b=B):
    """Get the point on the ellipsoid under the *query_point*."""
    x, y, z = query_point

    lat = geodetic_lat(query_point)
    lon = np.arctan2(y, x)
    e2_ = (a * a - b * b) / (a * a)
    n__ = a / np.sqrt(1 - e2_ * np.sin(lat)**2)
    nx_ = n__ * np.cos(lat) * np.cos(lon)
    ny_ = n__ * np.cos(lat) * np.sin(lon)
    nz_ = (1 - e2_) * n__ * np.sin(lat)

    return np.stack([nx_, ny_, nz_], axis=0)


def compute_yaw_steering(pos, vel):
    """Compute the yaw steering angle to compensate for Earth's rotation.

    Args:
        pos: Satellite position as column vector(s) in km (ECI frame).
        vel: Satellite velocity as column vector(s) in km/s (ECI frame).

    Returns:
        Yaw steering angle in radians.
    """
    r = vnorm(pos)
    v = vnorm(vel)
    lat = np.arcsin(pos[2] / r)
    return np.arctan2(OMEGA_EARTH * A * np.cos(lat), v)


def _local_frame(pos, vel):
    """Compute the satellite's local orbital reference frame.

    Returns (nadir, along_track, cross_track) as unit column vectors.
    """
    nadir = subpoint(-pos)
    nadir /= vnorm(nadir)
    along_track = vel / vnorm(vel)
    cross_track = np.cross(nadir, vel, 0, 0, 0)
    cross_track /= vnorm(cross_track)
    return nadir, along_track, cross_track


def _effective_yaw(yaw, yaw_steering, pos, vel, fovs_shape):
    """Compute the effective yaw angle, optionally adding the yaw-steering term."""
    if yaw_steering:
        yaw = yaw + compute_yaw_steering(pos, vel)
    if np.shape(yaw):
        yaw = np.broadcast_to(yaw, fovs_shape)
    return yaw


class ScanGeometry(object):
    """Description of the geometry of an instrument.

    *fovs* is the x and y viewing angles of the instrument. y is zero if the we
    talk about scanlines of course. *times* is the time of viewing of each
    angle relative to the start of the scanning, so it should have the same
    size as the *fovs*. *attitude* is the attitude correction to apply.
    """

    def __init__(self, fovs, times, attitude=(0, 0, 0)):
        """Initialize the class."""
        self.fovs = np.array(fovs)
        try:
            # assuming seconds
            self._times = np.asanyarray(times) * np.timedelta64(1000000000, "ns")
        except TypeError:
            self._times = np.asanyarray(times).astype("timedelta64[ns]")
        self.attitude = attitude

    def vectors(self, pos, vel, roll=0.0, pitch=0.0, yaw=0.0, yaw_steering=False):
        """Get unit vectors pointing to the different pixels.

        *pos* and *vel* are column vectors, or matrices of column
        vectors. Returns vectors as stacked rows.

        If *yaw_steering* is True, the yaw angle is computed from the
        satellite position and velocity to compensate for Earth's rotation.
        This is added to any explicit *yaw* value.
        """
        nadir, along_track, cross_track = _local_frame(pos, vel)
        effective_yaw = _effective_yaw(yaw, yaw_steering, pos, vel, self.fovs[0].shape)
        along_track_rotated = qrotate(nadir, along_track, self.fovs[0] + roll)
        both_rotated = qrotate(along_track_rotated, cross_track, self.fovs[1] + pitch)
        return qrotate(both_rotated, nadir, effective_yaw)

    def times(self, start_of_scan):
        """Return an array with the times of each scan line."""
        try:
            return np.array(self._times) + np.datetime64(start_of_scan)
        except ValueError:
            return np.array(self._times) + start_of_scan


class Quaternion(object):
    """Some class, that I don't know what is doing..."""

    def __init__(self, scalar, vector):
        """Initialize the class."""
        self.__x, self.__y, self.__z = vector.reshape((3, -1))
        self.__w = scalar.ravel()

    def rotation_matrix(self):
        """Get the rotation matrix."""
        x, y, z, w = self.__x, self.__y, self.__z, self.__w
        zero = np.zeros_like(x)
        return np.array(
            ((w**2 + x**2 - y**2 - z**2,
              2 * x * y + 2 * z * w,
              2 * x * z - 2 * y * w,
              zero),
             (2 * x * y - 2 * z * w,
              w**2 - x**2 + y**2 - z**2,
              2 * y * z + 2 * x * w,
              zero),
             (2 * x * z + 2 * y * w,
              2 * y * z - 2 * x * w,
              w**2 - x**2 - y**2 + z**2,
              zero),
             (zero, zero, zero, w**2 + x**2 + y**2 + z**2)))


def qrotate(vector, axis, angle):
    """Rotate *vector* around *axis* by *angle* (in radians).

    *vector* is a matrix of column vectors, as is *axis*.
    This function uses quaternion rotation.
    """
    n_axis = axis / vnorm(axis)
    sin_angle = np.expand_dims(np.sin(angle / 2), 0)
    if np.ndim(n_axis) == 1:
        n_axis = np.expand_dims(n_axis, 1)
        p__ = np.dot(n_axis, sin_angle)[:, np.newaxis]
    else:
        p__ = n_axis * sin_angle

    q__ = Quaternion(np.cos(angle / 2), p__)
    shape = vector.shape
    return np.einsum("kj, ikj->ij",
                     vector.reshape((3, -1)),
                     q__.rotation_matrix()[:3, :3]).reshape(shape)


def get_lonlatalt(pos, utc_time):
    """Calculate sublon, sublat and altitude of satellite, considering the earth an ellipsoid."""
    trf = Transformer.from_crs(dict(proj="geocent"), dict(proj="latlong"))
    lon, lat, alt = trf.transform(*(pos * 1000))
    lon = lon - np.rad2deg(astronomy.gmst(utc_time))
    lon = np.where(lon < -180, lon + 360, lon)
    return lon, lat, alt

# END OF DIRTY STUFF


def compute_pixels(orb, sgeom, times, rpy=(0.0, 0.0, 0.0), yaw_steering=False):
    """Compute cartesian coordinates of the pixels in instrument scan."""
    if isinstance(orb, (list, tuple)):
        tle1, tle2 = orb
        orb = Orbital("mysatellite", line1=tle1, line2=tle2)

    # get position and velocity for each time of each pixel
    pos, vel = orb.get_position(times, normalize=False)

    # now, get the vectors pointing to each pixel
    vectors = sgeom.vectors(pos, vel, *rpy, yaw_steering=yaw_steering)

    # compute intersection of lines (directed by vectors and passing through
    # (0, 0, 0)) and ellipsoid. Derived from:
    # http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

    # do the computation between line and ellipsoid (WGS 84)
    # NB: AAPP uses GRS 80...
    centre = -pos
    a__ = 6378.137  # km
    # b__ = 6356.75231414 # km, GRS80
    b__ = 6356.752314245  # km, WGS84
    radius = np.array([[1 / a__, 1 / a__, 1 / b__]]).T
    shape = vectors.shape

    xr_ = vectors.reshape([3, -1]) * radius
    cr_ = centre.reshape([3, -1]) * radius
    ldotc = np.einsum("ij,ij->j", xr_, cr_)
    lsq = np.einsum("ij,ij->j", xr_, xr_)
    csq = np.einsum("ij,ij->j", cr_, cr_)

    d1_ = (ldotc - np.sqrt(ldotc ** 2 - csq * lsq + lsq)) / lsq

    # return the actual pixel positions
    return vectors * d1_.reshape(shape[1:]) - centre


def norm(v):
    """Return the norm of the vector *v*."""
    return np.sqrt(np.dot(v, v.conj()))


def mnorm(m, axis=None):
    """Norm of a matrix of vectors stacked along the *axis* dimension."""
    if axis is None:
        axis = np.ndim(m) - 1
    return np.sqrt((m**2).sum(axis))


def vnorm(m):
    """Norms of a matrix of column vectors."""
    return np.sqrt((m**2).sum(0))


def hnorm(m):
    """Norms of a matrix of row vectors."""
    return np.sqrt((m**2).sum(1))


def _n_scans_from_duration(duration, time_sampling):
    """Compute the number of scan lines covering a given duration."""
    if time_sampling <= np.timedelta64(0):
        return 1
    return max(1, int(duration / time_sampling))


def _sample_indices(count, n_samples):
    """Return n_samples evenly-spaced integer indices from 0 to count-1."""
    return np.linspace(0, count - 1, n_samples, dtype=int)


def _boundary_scan_pixel_pairs(scan_indices, pixel_indices):
    """Build ordered (scan, pixel) pairs tracing the swath boundary polygon.

    Traverses: bottom edge left→right, right edge bottom→top,
    top edge right→left, left edge top→bottom, closing point.
    """
    scans, pixels = [], []

    for px in pixel_indices:               # bottom edge: left → right
        scans.append(scan_indices[0])
        pixels.append(px)
    for sc in scan_indices[1:]:            # right edge: bottom → top
        scans.append(sc)
        pixels.append(pixel_indices[-1])
    for px in reversed(pixel_indices[:-1]): # top edge: right → left
        scans.append(scan_indices[-1])
        pixels.append(px)
    for sc in reversed(scan_indices[1:-1]): # left edge: top → bottom
        scans.append(sc)
        pixels.append(pixel_indices[0])

    scans.append(scans[0])     # close the polygon
    pixels.append(pixels[0])
    return np.array(scans), np.array(pixels)


def _geolocate_boundary(swath, edge_scans, edge_pixels, start_time, tle, rpy):
    """Geolocate a set of (scan, pixel) boundary pairs and return (lons, lats)."""
    x_fovs, y_fovs = swath.scanline.angles(edge_pixels)
    sgeom = ScanGeometry(np.vstack((x_fovs, y_fovs)), edge_scans * swath.time_sampling)
    s_times = sgeom.times(start_time)
    pixels_pos = compute_pixels(tle, sgeom, s_times, rpy)
    lons, lats, _ = get_lonlatalt(pixels_pos, s_times)
    return lons, lats


def bounding_box(swath, start_time, end_time, tle, points_per_edge=10, rpy=(0.0, 0.0, 0.0)):
    """Compute a bounding polygon for a satellite swath.

    Args:
        swath: A PushbroomSwath with scanline and time_sampling attributes.
        start_time: Start of the observation period (datetime).
        end_time: End of the observation period (datetime).
        tle: TLE data as (line1, line2) tuple.
        points_per_edge: Number of sample points per edge (including corners).
        rpy: Roll, pitch, yaw corrections.

    Returns:
        Tuple of (lons, lats) arrays forming a closed polygon.
    """
    duration = np.datetime64(end_time) - np.datetime64(start_time)
    n_scans = _n_scans_from_duration(duration, swath.time_sampling)

    scan_indices = _sample_indices(n_scans, points_per_edge)
    pixel_indices = _sample_indices(swath.scanline.pixels_per_scan, points_per_edge)

    edge_scans, edge_pixels = _boundary_scan_pixel_pairs(scan_indices, pixel_indices)
    return _geolocate_boundary(swath, edge_scans, edge_pixels, start_time, tle, rpy)


if __name__ == "__main__":
    # NOAA 18 (from the 2011-10-12, 16:55 utc)
    # 1 28654U 05018A   11284.35271227  .00000478  00000-0  28778-3 0  9246
    # 2 28654  99.0096 235.8581 0014859 135.4286 224.8087 14.11526826329313

    noaa18_tle1 = "1 28654U 05018A   11284.35271227  .00000478  00000-0  28778-3 0  9246"
    noaa18_tle2 = "2 28654  99.0096 235.8581 0014859 135.4286 224.8087 14.11526826329313"

    import datetime as dt
    t = dt.datetime(2011, 10, 12, 13, 45)

    # edge and centre of an avhrr scanline
    # sgeom = ScanGeometry([(-0.9664123687741623, 0),
    #                      (0, 0)],
    #                     [0, 0.0, ])
    # print compute_pixels((noaa18_tle1, noaa18_tle2), sgeom, t)

    # avhrr swath
    scanline_nb = 1

    # building the avhrr angles, 2048 pixels from +55.37 to -55.37 degrees
    avhrr = np.vstack(((np.arange(2048) - 1023.5) / 1024 * np.deg2rad(-55.37),
                       np.zeros((2048,)))).transpose()
    avhrr = np.tile(avhrr, [scanline_nb, 1])
    # building the corresponding times array
    offset = np.arange(scanline_nb) * 0.1667
    times = (np.tile(np.arange(2048) * 0.000025 + 0.0025415, [scanline_nb, 1])
             + np.expand_dims(offset, 1))
    # build the scan geometry object
    sgeom = ScanGeometry(avhrr, times.ravel())

    # print the lonlats for the pixel positions
    s_times = sgeom.times(t)
    pixels_pos = compute_pixels((noaa18_tle1, noaa18_tle2), sgeom, s_times)
    print(get_lonlatalt(pixels_pos, s_times))
