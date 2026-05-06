"""Some instrument definitions to use with geoloc.

To define an instrument, one must first define the scan angles (in radians)
around x (along-track vector) and y (cross-track vector). the y scan angles are
just 0 in the case of scanline based instruments (like avhrr), but can be
different if the instrument is forward and/or backward scanning (e.g. viirs or
modis).

For the instrument to be defined completely, one must also provide the
observation times (in seconds, respective to the nominal scan time) for the
different pixels.

Both scan angles and scan times are then combined into a ScanGeometry object.
"""

import warnings
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from pyorbital.geoloc import ScanGeometry


@dataclass
class SingleLinePushbroomScan:
    """Definition of a single-line pushbroom instrument scan geometry.

    Args:
        left_angle: Scan angle at the left edge of the swath (degrees).
        right_angle: Scan angle at the right edge of the swath (degrees).
        pixels_per_scan: Number of pixels per scan line.
        forward_angle: Along-track viewing angle (degrees), default 0.
    """

    left_angle: float
    right_angle: float
    pixels_per_scan: int
    forward_angle: float = 0

    def angles(self, pixels: slice | ArrayLike | None = None):
        """Compute the scan angles for the given pixel positions.

        Args:
            pixels: Pixel positions to compute angles for. Can be a slice,
                an array-like of indices, or None for all pixels.

        Returns:
            Tuple of (cross_track_angles, along_track_angles) in radians.
        """
        left, right, count = self._resolve_pixel_range(pixels)
        x_fovs = np.linspace(np.deg2rad(left), np.deg2rad(right), count)
        y_fovs = np.full(count, np.deg2rad(self.forward_angle))
        return x_fovs, y_fovs

    def _resolve_pixel_range(self, pixels):
        """Resolve a pixel selector to (left_angle_deg, right_angle_deg, count)."""
        if pixels is None:
            return self.left_angle, self.right_angle, self.pixels_per_scan
        positions = self._pixel_positions(pixels)
        scale = (self.right_angle - self.left_angle) / (self.pixels_per_scan - 1)
        left = positions[0] * scale + self.left_angle
        right = positions[-1] * scale + self.left_angle
        return left, right, len(positions)

    def _pixel_positions(self, pixels):
        """Convert a pixel selector (slice or array-like) to a sequence of positions."""
        try:
            return range(*pixels.indices(self.pixels_per_scan))
        except AttributeError:
            return pixels


@dataclass
class PushbroomSwath:
    """A pushbroom swath combining a scanline definition with time sampling.

    Args:
        scanline: The instrument scanline definition.
        time_sampling: Time interval between consecutive scan lines.
    """

    scanline: SingleLinePushbroomScan
    time_sampling: np.timedelta64

    def scan_geometry(self, scan_lines: slice, pixels: slice | ArrayLike | None = None):
        """Generate a ScanGeometry for the given scan lines and pixel subset.

        Args:
            scan_lines: Slice selecting which scan line numbers to include.
                The slice stop is required; start defaults to 0, step to 1.
            pixels: Optional pixel subset (slice or array of indices).

        Returns:
            A ScanGeometry object with the appropriate fovs and times.
        """
        scans = range(*scan_lines.indices(scan_lines.stop))
        x_fovs, y_fovs = self.scanline.angles(pixels)
        fovs = self._tiled_fovs(x_fovs, y_fovs, len(scans))
        times = self._scan_times(scans, len(x_fovs))
        return ScanGeometry(fovs, times)

    def _tiled_fovs(self, x_fovs, y_fovs, n_scans):
        """Stack and tile scan angles across all scan lines."""
        one_line = np.vstack((x_fovs, y_fovs))
        return np.tile(one_line[:, np.newaxis, :], [1, n_scans, 1])

    def _scan_times(self, scans, pixels_len):
        """Build a (n_scans × pixels_len) time offset array."""
        times = np.zeros((len(scans), pixels_len), dtype=self.time_sampling.dtype)
        times += np.array(scans)[:, np.newaxis] * self.time_sampling
        return times


@dataclass
class WhiskbroomScan:
    """Definition of a cross-track (whiskbroom) scanning instrument.

    Pixels are acquired sequentially across track, one line at a time.

    Args:
        pixels_per_scan: Number of pixels per scan line.
        scan_angle: Half-swath angle in degrees (positive). Pixels run
            from +scan_angle (left edge) to -scan_angle (right edge).
        scan_rate: Duration of one complete scan cycle in seconds.
        pixel_dwell_time: Time per pixel measurement in seconds.
        sync_time: Delay before the first pixel measurement in seconds.
    """

    pixels_per_scan: int
    scan_angle: float
    scan_rate: float
    pixel_dwell_time: float
    sync_time: float = 0.0

    def angles(self, scan_points=None):
        """Compute cross-track angles for the given pixel positions.

        Args:
            scan_points: Pixel indices (array-like or None for all pixels).

        Returns:
            Cross-track angles in radians, running from +scan_angle to
            -scan_angle as pixel index increases.
        """
        if scan_points is None:
            scan_points = np.arange(self.pixels_per_scan)
        scan_points = np.asanyarray(scan_points)
        return (scan_points / (self.pixels_per_scan * 0.5 - 0.5) - 1) * np.deg2rad(-self.scan_angle)

    def scan_geometry(self, scans_nb, scan_points=None):
        """Generate a ScanGeometry for the given number of scan lines.

        Args:
            scans_nb: Number of scan lines.
            scan_points: Optional subset of pixel indices.

        Returns:
            A ScanGeometry object with the appropriate fovs and times.
        """
        scan_points = self._resolve_scan_points(scan_points)
        samples = self._tiled_samples(scan_points, int(scans_nb))
        times = self._scan_times(scan_points, int(scans_nb))
        return ScanGeometry(samples, times)

    def _resolve_scan_points(self, scan_points):
        """Return scan_points as a numpy array, defaulting to all pixels."""
        if scan_points is None:
            return np.arange(self.pixels_per_scan)
        return np.asanyarray(scan_points)

    def _tiled_samples(self, scan_points, n_scans):
        """Stack cross-track and along-track angles, tiled across scan lines."""
        cross_track = self.angles(scan_points)
        one_line = np.vstack((cross_track, np.zeros(len(scan_points))))
        return np.tile(one_line[:, np.newaxis, :], [1, n_scans, 1])

    def _scan_times(self, scan_points, n_scans):
        """Build (n_scans × n_pixels) time array with per-pixel and per-scan offsets."""
        pixel_times = scan_points * self.pixel_dwell_time + self.sync_time
        scan_offsets = np.arange(n_scans) * self.scan_rate
        return np.tile(pixel_times, [n_scans, 1]) + np.expand_dims(scan_offsets, 1)


AMSU_A_SCAN = WhiskbroomScan(
    pixels_per_scan=30,
    scan_angle=48.3,
    scan_rate=8.0,
    pixel_dwell_time=0.2,
    sync_time=0.00355,
)

MHS_SCAN = WhiskbroomScan(
    pixels_per_scan=90,
    scan_angle=49.444,
    scan_rate=8 / 3.0,
    pixel_dwell_time=(8 / 3.0 - 1) / 90.0,
    sync_time=0.0,
)

HIRS4_SCAN = WhiskbroomScan(
    pixels_per_scan=56,
    scan_angle=49.5,
    scan_rate=6.4,
    pixel_dwell_time=6.4 / 56.0,
    sync_time=0.0,
)

ATMS_SCAN = WhiskbroomScan(
    pixels_per_scan=96,
    scan_angle=52.7,
    scan_rate=8 / 3.0,
    pixel_dwell_time=18e-3,
    sync_time=0.0,
)

MWHS2_SCAN = WhiskbroomScan(
    pixels_per_scan=98,
    scan_angle=53.35,
    scan_rate=8 / 3.0,
    pixel_dwell_time=(8 / 3.0 - 1) / 98.0,
    sync_time=0.0,
)


@dataclass
class MultiLineWhiskbroomScan:
    """Definition of a multi-line cross-track (whiskbroom) scanning instrument.

    Each sweep of the cross-track mirror captures *lines_per_scan* simultaneous
    lines via a detector array stacked along the flight direction.  Examples:
    MERSI-2, MERSI-3, MODIS, VIIRS (all have multiple detector rows per scan).

    Args:
        pixels_per_scan: Number of pixels per scan line (across-track).
        scan_angle: Half-swath angle in degrees (positive). Pixels run
            from +scan_angle (left edge) to -scan_angle (right edge).
        scan_rate: Duration of one complete scan cycle in seconds.
        pixel_dwell_time: Time per pixel measurement in seconds.
        lines_per_scan: Number of detector rows along-track captured per sweep.
        along_track_step: Angular spacing between detector rows in radians.
            Typically pixel_size_km / orbit_altitude_km.
        sync_time: Delay before the first pixel measurement in seconds.
    """

    pixels_per_scan: int
    scan_angle: float
    scan_rate: float
    pixel_dwell_time: float
    lines_per_scan: int
    along_track_step: float
    sync_time: float = 0.0

    def cross_track_angles(self, scan_points=None):
        """Compute cross-track angles for the given pixel positions.

        Returns:
            Cross-track angles in radians, running from +scan_angle to
            -scan_angle as pixel index increases.
        """
        if scan_points is None:
            scan_points = np.arange(self.pixels_per_scan)
        scan_points = np.asanyarray(scan_points)
        return (scan_points / (self.pixels_per_scan * 0.5 - 0.5) - 1) * np.deg2rad(-self.scan_angle)

    def along_track_angles(self):
        """Compute along-track angles for each detector row.

        Returns:
            Along-track angles in radians, symmetric around zero, one per row.
        """
        rows = np.arange(self.lines_per_scan)
        return ((self.lines_per_scan - 1) / 2.0 - rows) * self.along_track_step

    def scan_geometry(self, n_scans, scan_points=None):
        """Generate a ScanGeometry for the given number of complete scans.

        Each scan produces *lines_per_scan* output lines.

        Args:
            n_scans: Number of scan cycles.
            scan_points: Optional subset of pixel indices.

        Returns:
            A ScanGeometry with fovs shape ``(2, n_scans*lines_per_scan, n_pixels)``
            and times shape ``(n_scans*lines_per_scan, n_pixels)``.
        """
        if scan_points is None:
            scan_points = np.arange(self.pixels_per_scan)
        scan_points = np.asanyarray(scan_points)
        n_scans = int(n_scans)
        fovs = self._build_fovs(scan_points, n_scans)
        times = self._build_times(scan_points, n_scans)
        return ScanGeometry(fovs, times, lines_per_scan=self.lines_per_scan)

    def _build_fovs(self, scan_points, n_scans):
        """Build (2, n_scans*L, n_pixels) FOV array."""
        L = self.lines_per_scan
        cross = self.cross_track_angles(scan_points)           # (N,)
        along = self.along_track_angles()                      # (L,)
        # cross-track: same for all L lines within each scan
        cross_tiled = np.tile(cross, (n_scans * L, 1))        # (n_scans*L, N)
        # along-track: repeating L-row block, constant across pixels
        along_block = np.tile(along[:, np.newaxis], (1, len(scan_points)))  # (L, N)
        along_tiled = np.tile(along_block, (n_scans, 1))      # (n_scans*L, N)
        return np.stack([cross_tiled, along_tiled])            # (2, n_scans*L, N)

    def _build_times(self, scan_points, n_scans):
        """Build (n_scans*L, n_pixels) time array.

        All detector rows within the same scan share identical per-pixel times
        (the mirror position is the same for all rows at each instant).
        """
        L = self.lines_per_scan
        pixel_times = scan_points * self.pixel_dwell_time + self.sync_time  # (N,)
        scan_offsets = np.arange(n_scans) * self.scan_rate                  # (n_scans,)
        # one row per scan line (n_scans*L rows total, L identical per scan)
        per_scan_times = (np.tile(pixel_times, (n_scans, 1))
                          + scan_offsets[:, np.newaxis])                    # (n_scans, N)
        return np.repeat(per_scan_times, L, axis=0)                        # (n_scans*L, N)


# FY-3 MERSI-2 / MERSI-3 (same scan geometry for both series)
# Sources: WMO OSCAR; scan rate from EV_start_time in L1B data (1.5 s/scan)
# Calibrated against FY-3F NSMC L1B reference data (2026-03-20 08:25 pass):
#   scan_angle  = 55.1349° (WMO OSCAR lists 55.4° which is ~0.27° off)
#   along_track = 1/883 rad (orbital altitude 883 km, not 850)
#   sync_time   = −(381 pixels × 1.48 s / 2048) = −0.2753 s
#     The raw VCDU timestamp marks pixel 381 of the EV sweep (not pixel 0).
_MERSI_SWEEP_TIME = 1.48  # seconds of active sweep per scan
_MERSI_SYNC_TIME = -(381 * _MERSI_SWEEP_TIME / 2048)   # ≈ −0.2753 s
_MERSI_ALTITUDE_KM = 830.0                              # calibrated orbital altitude

MERSI_1KM_SCAN = MultiLineWhiskbroomScan(
    pixels_per_scan=2048,
    scan_angle=55.1349,
    scan_rate=1.5,
    pixel_dwell_time=_MERSI_SWEEP_TIME / 2048,
    lines_per_scan=10,
    along_track_step= 1 / _MERSI_ALTITUDE_KM,
    sync_time=_MERSI_SYNC_TIME,
)

MERSI_250M_SCAN = MultiLineWhiskbroomScan(
    pixels_per_scan=8192,
    scan_angle=55.1349,
    scan_rate=1.5,
    pixel_dwell_time=_MERSI_SWEEP_TIME / 8192,
    lines_per_scan=40,
    along_track_step=0.25 / _MERSI_ALTITUDE_KM,
    sync_time=_MERSI_SYNC_TIME,
)


@dataclass
class FocalPlaneWhiskbroomScan:
    """Multi-line whiskbroom scan with explicit focal-plane geometry and boresight.

    Extends the capability of :class:`MultiLineWhiskbroomScan` by deriving
    per-detector along-track angles directly from the instrument's focal-plane
    geometry, and by applying an optional boresight rotation matrix that maps
    the instrument frame to the spacecraft body frame.

    This is the appropriate model for instruments like FY-3 MERSI where
    different spectral bands sit at different along-track positions in the focal
    plane, causing per-band geolocation differences that the simpler
    ``along_track_step`` approximation cannot capture.

    The focal-plane model follows the reference implementation in
    CalTelescoppViewVector1KM (mepgps_F3D.c).  For detector row *i* (0-based,
    forward-looking rows first)::

        y_i = ((lines_per_scan - 1) / 2 - i) * det_space + det_position[1]
        x_i = det_position[0]
        z_i = focal_length
        unit_vec = [x_i, y_i, z_i] / norm([x_i, y_i, z_i])

    All length units (``focal_length``, ``det_space``, ``det_position``) must
    be consistent (e.g. all mm); only the ratios matter because vectors are
    normalised.

    Args:
        pixels_per_scan: Number of pixels per scan line (across-track).
        scan_angle: Half-swath angle in degrees (positive). Pixels run from
            +scan_angle (left edge) to −scan_angle (right edge).
        scan_rate: Duration of one complete scan cycle in seconds.
        pixel_dwell_time: Time per pixel measurement in seconds.
        lines_per_scan: Number of detector rows along-track per sweep.
        focal_length: Telescope focal length in consistent length units.
        det_space: Along-track detector pitch in the same length units.
        sync_time: Delay before the first pixel measurement in seconds.
        det_position: ``(x, y)`` band-centre offset in the focal plane,
            in the same length units.  *x* is the cross-track offset,
            *y* is the along-track offset.  Defaults to ``(0.0, 0.0)``.
        boresight: Optional 3×3 rotation matrix **T_inst2sc** that maps the
            instrument telescope frame to the spacecraft body frame.  When
            ``None`` the identity is used.
    """

    pixels_per_scan: int
    scan_angle: float
    scan_rate: float
    pixel_dwell_time: float
    lines_per_scan: int
    focal_length: float
    det_space: float
    sync_time: float = 0.0
    det_position: tuple = (0.0, 0.0)
    boresight: "np.ndarray | None" = None

    def _focal_plane_unit_vectors(self):
        """Return (3, L) unit vectors in the instrument telescope frame."""
        rows = np.arange(self.lines_per_scan, dtype=float)
        y = (self.lines_per_scan - 1) / 2.0 - rows
        y = y * self.det_space + self.det_position[1]
        x = np.full(self.lines_per_scan, self.det_position[0], dtype=float)
        z = np.full(self.lines_per_scan, self.focal_length, dtype=float)
        vecs = np.stack([x, y, z])                         # (3, L)
        return vecs / np.linalg.norm(vecs, axis=0)         # normalised

    def _body_frame_vectors(self):
        """Return (3, L) unit vectors in the spacecraft body frame."""
        vecs = self._focal_plane_unit_vectors()
        if self.boresight is not None:
            vecs = np.asarray(self.boresight, dtype=float) @ vecs
        return vecs

    def along_track_angles(self):
        """Compute along-track angles for each detector row.

        Derived from the focal-plane geometry and the boresight matrix.
        Positive angles point in the forward flight direction.

        Returns:
            Along-track angles in radians, one per detector row.
        """
        vecs = self._body_frame_vectors()                  # (3, L)
        # Convention: axis 1 = along-track, axis 2 = nadir/boresight
        return np.arctan2(vecs[1], vecs[2])

    def cross_track_angles(self, scan_points=None):
        """Compute cross-track angles for the given pixel positions.

        The standard linear scan-angle formula is used, with a scalar
        band-centre cross-track offset derived from the boresighted focal-plane
        geometry (``det_position[0]`` after applying the boresight matrix).

        Returns:
            Cross-track angles in radians, running from +scan_angle to
            −scan_angle as pixel index increases.
        """
        if scan_points is None:
            scan_points = np.arange(self.pixels_per_scan)
        scan_points = np.asanyarray(scan_points)
        base = (scan_points / (self.pixels_per_scan * 0.5 - 0.5) - 1) * np.deg2rad(-self.scan_angle)
        vecs = self._body_frame_vectors()                  # (3, L)
        # Mean cross-track offset of the band centre across all detector rows
        ct_offset = float(np.arctan2(vecs[0], vecs[2]).mean())
        return base + ct_offset

    def scan_geometry(self, n_scans, scan_points=None):
        """Generate a ScanGeometry for the given number of complete scans.

        Each scan produces *lines_per_scan* output lines.

        Args:
            n_scans: Number of scan cycles.
            scan_points: Optional subset of pixel indices.

        Returns:
            A :class:`~pyorbital.geoloc.ScanGeometry` with ``fovs`` shape
            ``(2, n_scans*lines_per_scan, n_pixels)`` and ``times`` shape
            ``(n_scans*lines_per_scan, n_pixels)``.
        """
        if scan_points is None:
            scan_points = np.arange(self.pixels_per_scan)
        scan_points = np.asanyarray(scan_points)
        n_scans = int(n_scans)
        fovs = self._build_fovs(scan_points, n_scans)
        times = self._build_times(scan_points, n_scans)
        return ScanGeometry(fovs, times, lines_per_scan=self.lines_per_scan)

    def _build_fovs(self, scan_points, n_scans):
        """Build (2, n_scans*L, n_pixels) FOV array."""
        L = self.lines_per_scan
        cross = self.cross_track_angles(scan_points)            # (N,)
        along = self.along_track_angles()                       # (L,)
        cross_tiled = np.tile(cross, (n_scans * L, 1))         # (n_scans*L, N)
        along_block = np.tile(along[:, np.newaxis], (1, len(scan_points)))  # (L, N)
        along_tiled = np.tile(along_block, (n_scans, 1))       # (n_scans*L, N)
        return np.stack([cross_tiled, along_tiled])             # (2, n_scans*L, N)

    def _build_times(self, scan_points, n_scans):
        """Build (n_scans*L, n_pixels) time array.

        All detector rows within the same scan share identical per-pixel times
        (the mirror position is the same for all rows at each instant).
        """
        L = self.lines_per_scan
        pixel_times = scan_points * self.pixel_dwell_time + self.sync_time  # (N,)
        scan_offsets = np.arange(n_scans) * self.scan_rate                  # (n_scans,)
        per_scan_times = (np.tile(pixel_times, (n_scans, 1))
                          + scan_offsets[:, np.newaxis])                    # (n_scans, N)
        return np.repeat(per_scan_times, L, axis=0)                         # (n_scans*L, N)


################################################################
#
#   AVHRR
#
################################################################


def avhrr(scans_nb, scan_points,
          scan_angle=55.37, frequency=1 / 6.0, apply_offset=True):
    """Definition of the avhrr instrument.

    Source: NOAA KLM User's Guide, Appendix J
    http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/j/app-j.htm
    """
    # build the avhrr instrument (scan angles)
    avhrr_inst = np.vstack(((scan_points / 1023.5 - 1)
                            * np.deg2rad(-scan_angle),
                            np.zeros((len(scan_points),))))

    avhrr_inst = np.tile(
        avhrr_inst[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

    # building the corresponding times array
    # times = (np.tile(scan_points * 0.000025 + 0.0025415, [scans_nb, 1])
    #         + np.expand_dims(offset, 1))

    times = np.tile(scan_points * 0.000025, [np.int32(scans_nb), 1])
    if apply_offset:
        offset = np.arange(np.int32(scans_nb)) * frequency
        times += np.expand_dims(offset, 1)

    return ScanGeometry(avhrr_inst, times)


def avhrr_gac(scan_times, scan_points,
              scan_angle=55.37, frequency=0.5):
    """Definition of the avhrr instrument, gac version.

    Source: NOAA KLM User's Guide, Appendix J
    http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/j/app-j.htm
    """
    warnings.warn("avhrr_gac is replaced with avhrr_from_times or avhrr_gac_from_times",
                  DeprecationWarning)
    try:
        offset = np.array([(t - scan_times[0]).seconds +
                           (t - scan_times[0]).microseconds / 1000000.0 for t in scan_times])
    except TypeError:
        offset = np.arange(scan_times) * frequency
    scans_nb = len(offset)

    avhrr_inst = np.vstack(((scan_points / 1023.5 - 1)
                            * np.deg2rad(-scan_angle),
                            np.zeros((len(scan_points),))))

    avhrr_inst = np.tile(
        avhrr_inst[:, np.newaxis, :], [1, np.int32(scans_nb), 1])
    # building the corresponding times array
    times = (np.tile(scan_points * 0.000025, [scans_nb, 1])
             + np.expand_dims(offset, 1))
    return ScanGeometry(avhrr_inst, times)


def avhrr_from_times(scan_times, scan_points, scan_angle=55.37):
    """Definition of the avhrr instrument.

    Source: NOAA KLM User's Guide, Appendix J
    http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/j/app-j.htm

    :scan_times: Observation times (datetime object)
    :scan_points: Across track pixel positions
    :scan_angle: Maximum scan angle of the outermost FOV
    """
    offset = np.array([(t - scan_times[0]).total_seconds() for t in scan_times])
    scan_points = np.asanyarray(scan_points)
    scans_nb = len(offset)

    avhrr_inst = np.vstack([(scan_points / 1023.5 - 1) * np.deg2rad(-scan_angle),
                            np.zeros(len(scan_points))])

    avhrr_inst = np.tile(avhrr_inst[:, np.newaxis, :], [1, scans_nb, 1])
    # building the corresponding times array
    times = (np.tile(scan_points * 0.000025, [scans_nb, 1])
             + np.expand_dims(offset, 1))
    return ScanGeometry(avhrr_inst, times)


def avhrr_gac_from_times(scan_times, scan_points, scan_angle=55.37):
    """Definition of the avhrr instrument, gac version.

    Source: NOAA KLM User's Guide, Appendix J
    http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/j/app-j.htm

    :scan_times: Observation times (datetime object)
    :scan_points: Across track pixel positions
    :scan_angle: Maximum scan angle of the outermost FOV
    """
    offset = np.array([(t - scan_times[0]).total_seconds() for t in scan_times])
    scan_points = np.asanyarray(scan_points)
    scans_nb = len(offset)

    # The AVHRR swath is +/- 55.4 degrees, which puts the outermost FOV at
    #   55.4 * 1023.5 / 1024 = 55.373
    # GAC pixels are the average of four FOVs with sampling: ..1111.2222.----.NNNN..
    # so we need to reduce the scan_angle by a factor (1023.5 - 3.5) / 1023.5
    # to calculate the GAC pixel centres
    gac_angle = scan_angle * (1023.5 - 3.5) / 1023.5

    avhrr_inst = np.vstack([(scan_points / 204 - 1) * np.deg2rad(-gac_angle),
                            np.zeros(len(scan_points))])

    avhrr_inst = np.tile(
        avhrr_inst[:, np.newaxis, :], [1, np.int32(scans_nb), 1])
    # building the corresponding times array
    times = (np.tile(scan_points * 0.000125, [scans_nb, 1])
             + np.expand_dims(offset, 1))
    return ScanGeometry(avhrr_inst, times)


################################################################
# avhrr, all pixels

# build the scan geometry object


def avhrr_all_geom(scans_nb):
    """Get all the AVHRR scan points."""
    # we take all pixels
    scan_points = np.arange(2048)
    return avhrr(scans_nb, scan_points)

################################################################
# avhrr, edge pixels

# build the scan geometry object


def avhrr_edge_geom(scans_nb):
    """Getting the AVHRR scan edges only."""
    scan_points = np.array([0, 2047])
    return avhrr(scans_nb, scan_points)

################################################################
# avhrr, every 40th pixel from the 24th (aapp style)

# build the scan geometry object


def avhrr_40_geom(scans_nb):
    """Description of the AVHRR scan in terms of every 40th pixel per line."""
    scan_points = np.arange(24, 2048, 40)
    return avhrr(scans_nb, scan_points)

################################################################
#
#   VIIRS
#
################################################################


def viirs(scans_nb, scan_indices=slice(0, None),
          chn_pixels=6400, scan_lines=32, scan_step=1):
    """Describe VIIRS instrument geometry, I-band by default.

    VIIRS scans several lines simultaneously (there are 16 detectors for each
    M-band, 32 detectors for each I-band) so the scan angles (and times) are
    two-dimensional arrays, contrary to AVHRR for example.

    scan_step: The increment in number of scans. E.g. if scan_step is 100 and
               the number of scans (scans_nb) is 10 then these 10 scans are
               distributed over the swath so that between each scan there are
               99 emtpy (excluded) scans

    """
    entire_width = np.arange(chn_pixels)
    scan_points = entire_width[scan_indices].astype("int")
    scan_pixels = len(scan_points)

    # Initial angle 55.84 deg replaced with 56.28 deg found in
    # VIIRS User's Guide from NESDIS, version 1.2 (09/10/2013).
    # Ref : NOAA Technical Report NESDIS 142.
    # Seems to be better (not quantified).
    across_track = \
        (scan_points / (chn_pixels / 2. - 0.5) - 1) * np.deg2rad(-56.28)
    y_max_angle = np.arctan2(11.87 / 2, 824.0)
    along_track = \
        -(np.arange(scan_lines) / (scan_lines / 2. - 0.5) - 1) * \
        y_max_angle
    scan = np.dstack((np.tile(across_track, (scan_lines, 1)).T,
                      np.tile(along_track, (scan_pixels, 1))))
    npp = np.tile(scan, [scans_nb, 1]).T

    # from the timestamp in the filenames, a granule takes 1:25.400 to record
    # (85.4 seconds) so 1.779166667 would be the duration of 1 scanline (48
    # scans per granule) dividing the duration of a single scan by a width of
    # 6400 pixels results in 0.0002779947917 seconds for each column of 32
    # pixels in the scanline

    # the individual times per pixel are probably wrong, unless the scanning
    # behaves the same as for AVHRR, The VIIRS sensor rotates to allow internal
    # calibration before each scanline. This would imply that the scanline
    # always moves in the same direction.  more info @
    # http://www.eoportal.org/directory/pres_NPOESSNationalPolarorbitingOperationalEnvironmentalSatelliteSystem.html

    SEC_EACH_SCANCOLUMN = 0.0002779947917
    sec_scan_duration = 1.779166667
    times = np.tile(scan_points * SEC_EACH_SCANCOLUMN,
                    [np.int32(scans_nb*scan_lines), 1])
    offset = np.repeat(np.arange(scans_nb) *
                       sec_scan_duration*scan_step, scan_lines)
    times += np.expand_dims(offset, 1)

    # build the scan geometry object
    return ScanGeometry(npp, times)


def viirs_edge_geom(scans_nb):
    """Definition of the VIIRS scane edges."""
    scan_indices = [0, -1]
    return viirs(scans_nb, scan_indices)


################################################################
#
#   AMSU-A
#
################################################################

def amsua(scans_nb, scan_points=None):
    """Describe AMSU-A instrument geometry.

    Parameters:
       scans_nb | int -  number of scan lines

     Keywords:
     * scan_points - FIXME!

    Returns:
       pyorbital.geoloc.ScanGeometry object

    """
    return AMSU_A_SCAN.scan_geometry(scans_nb, scan_points)


################################################################
#
#   MHS
#
################################################################

def mhs(scans_nb, scan_points=None):
    """Describe MHS instrument geometry.

    See:

    - https://www.eumetsat.int/website/home/Satellites/CurrentSatellites/Metop/MetopDesign/MHS/index.html
    - https://www1.ncdc.noaa.gov/pub/data/satellite/publications/podguides/
          N-15%20thru%20N-19/pdf/0.0%20NOAA%20KLM%20Users%20Guide.pdf
      (NOAA KLM Users Guide –August 2014 Revision)

    Parameters:
       scans_nb | int -  number of scan lines

     Keywords:
     * scan_points - FIXME!

    Returns:
       pyorbital.geoloc.ScanGeometry object

    """
    return MHS_SCAN.scan_geometry(scans_nb, scan_points)


################################################################
#
#   HIRS/4
#
################################################################

def hirs4(scans_nb, scan_points=None):
    """Describe HIRS/4 instrument geometry.

    See:
    - https://www.eumetsat.int/website/home/Satellites/CurrentSatellites/Metop/MetopDesign/HIRS/index.html
    - https://www1.ncdc.noaa.gov/pub/data/satellite/publications/podguides/
          N-15%20thru%20N-19/pdf/0.0%20NOAA%20KLM%20Users%20Guide.pdf
      (NOAA KLM Users Guide –August 2014 Revision)

    Parameters:
       scans_nb | int -  number of scan lines

     Keywords:
     * scan_points - FIXME!

    Returns:
       pyorbital.geoloc.ScanGeometry object

    """
    return HIRS4_SCAN.scan_geometry(scans_nb, scan_points)


################################################################
#
#   ATMS
#
################################################################

def atms(scans_nb, scan_points=None):
    """Describe ATMS instrument geometry.

    See:

    - https://dtcenter.org/com-GSI/users/docs/presentations/2013_workshop/
          Garrett_GSI_2013.pdf (Assimilation of Suomi-NPP ATMS, Kevin Garrett et al., August 8, 2013)
    - https://www.star.nesdis.noaa.gov/star/documents/meetings/2016JPSSAnnual/
          S4/S4_13_JPSSScience2016_session4Part2_ATMS_Scan_Reversal_HYANG.pdf
          (Suomi NPP ATMS Scan Reversal Study, Hu (Tiger) Yang, NOAA/STAR ATMS SDR Working Group)

    Parameters:
       scans_nb | int -  number of scan lines

     Keywords:
     * scan_points - FIXME!

    Returns:
       pyorbital.geoloc.ScanGeometry object

    """
    return ATMS_SCAN.scan_geometry(scans_nb, scan_points)


################################################################
#
#   MWHS-2
#
################################################################

def mwhs2(scans_nb, scan_points=None):
    """Describe MWHS-2 instrument geometry.

    The scanning period is 2.667 s. Main beams of the antenna scan over the ob-
    serving swath (±53.35◦ from nadir) in the cross-track direction at a
    constant time of 1.71 s. There are 98 pixels sampled per scan during 1.71s,
    and each sample has the same integration period.

    See:

       http://english.nssc.cas.cn/rh/rp/201501/W020150122580098790190.pdf

    Parameters:
       scans_nb | int -  number of scan lines

     Keywords:
     * scan_points - FIXME!

    Returns:
       pyorbital.geoloc.ScanGeometry object

    """
    return MWHS2_SCAN.scan_geometry(scans_nb, scan_points)


################################################################
#
#   OLCI
#
################################################################


OLCI_SCAN = SingleLinePushbroomScan(
    left_angle=46.5, right_angle=-22.1, pixels_per_scan=4000)

OLCI_SWATH = PushbroomSwath(
    scanline=OLCI_SCAN, time_sampling=np.timedelta64(44, "ms"))


def olci(scans_nb, scan_points=None, apply_offset=True):
    """Definition of the OLCI instrument.

    Source: Sentinel-3 OLCI Coverage
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-olci/coverage
    """
    if scan_points is not None:
        scan = SingleLinePushbroomScan(
            left_angle=OLCI_SCAN.left_angle,
            right_angle=OLCI_SCAN.right_angle,
            pixels_per_scan=len(scan_points))
    else:
        scan = OLCI_SCAN
    swath = PushbroomSwath(scanline=scan, time_sampling=OLCI_SWATH.time_sampling)
    geom = swath.scan_geometry(scan_lines=slice(int(scans_nb)))
    if not apply_offset:
        geom._times[:] = 0
    return geom


def ascat(scan_nb, scan_points=None):
    """Describing the ASCAT scanning geometry.

    make two scans one to the left and one to the right of the sub-satellite
    track.
    """
    if scan_points is None:
        scan_len = 42  # samples per scan
        scan_points = np.arange(42)
    else:
        scan_len = len(scan_points)

    scan_angle_inner = -25.0  # swath, degrees
    scan_angle_outer = -53.0  # swath, degrees
    scan_rate = 3.74747474747  # single scan, seconds
    if scan_len < 2:
        raise ValueError("Need at least two scan points!")

    sampling_interval = scan_rate / float(np.max(scan_points) + 1)

    # build the Metop/ascat instrument scan line angles
    scanline_angles_one = np.linspace(-np.deg2rad(scan_angle_outer),
                                      -np.deg2rad(scan_angle_inner), 21)
    scanline_angles_two = np.linspace(np.deg2rad(scan_angle_inner),
                                      np.deg2rad(scan_angle_outer), 21)

    scan_angles = np.concatenate(
        [scanline_angles_one, scanline_angles_two])[scan_points]

    inst = np.vstack((scan_angles, np.zeros(scan_len * 1,)))
    inst = np.tile(inst[:, np.newaxis, :], [1, np.int32(scan_nb), 1])

    # building the corresponding times array
    offset = np.arange(scan_nb) * scan_rate

    times = (np.tile(scan_points * sampling_interval,
                     [np.int32(scan_nb), 1]) + np.expand_dims(offset, 1))

    return ScanGeometry(inst, times)


SLSTR_NADIR_SCAN = SingleLinePushbroomScan(
    left_angle=46.5, right_angle=-22.1, pixels_per_scan=3000)


def slstr_nadir(scans_nb, scan_points=None):
    """Definition of the SLSTR instrument nadir swath.

    Source: Sentinel-3 SLSTR Coverage
    https://sentinel.esa.int/web/sentinel/user-guides/Sentinel-3-slstr/coverage
    """
    if scan_points is not None:
        scan = SingleLinePushbroomScan(
            left_angle=SLSTR_NADIR_SCAN.left_angle,
            right_angle=SLSTR_NADIR_SCAN.right_angle,
            pixels_per_scan=len(scan_points))
    else:
        scan = SLSTR_NADIR_SCAN
    swath = PushbroomSwath(scanline=scan, time_sampling=np.timedelta64(0))
    return swath.scan_geometry(scan_lines=slice(int(scans_nb)))
