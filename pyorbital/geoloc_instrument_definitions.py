#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013 - 2021 PyTroll Community

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Mikhail Itkin <itkin.m@gmail.com>
#   Adam Dybbroe <adam.dybbroe@smhi.se>

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

import numpy as np

from pyorbital.geoloc import ScanGeometry


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
    """Definition of the avhrr instrument, gac version

    Source: NOAA KLM User's Guide, Appendix J
    http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/j/app-j.htm
    """
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

################################################################
# avhrr, all pixels

# build the scan geometry object


def avhrr_all_geom(scans_nb):
    # we take all pixels
    scan_points = np.arange(2048)
    return avhrr(scans_nb, scan_points)

################################################################
# avhrr, edge pixels

# build the scan geometry object


def avhrr_edge_geom(scans_nb):
    # we take only edge pixels
    scan_points = np.array([0, 2047])
    return avhrr(scans_nb, scan_points)

################################################################
# avhrr, every 40th pixel from the 24th (aapp style)

# build the scan geometry object


def avhrr_40_geom(scans_nb):
    # we take only every 40th pixel
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
    scan_points = entire_width[scan_indices.astype('int')]
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
    # we take only edge pixels
    scan_indices = [0, -1]
    return viirs(scans_nb, scan_indices)


################################################################
#
#   AMSU-A
#
################################################################

def amsua(scans_nb, scan_points=None):
    """ Describe AMSU-A instrument geometry

    Parameters:
       scans_nb | int -  number of scan lines

     Keywords:
     * scan_points - FIXME!

    Returns:
       pyorbital.geoloc.ScanGeometry object

    """

    scan_len = 30  # 30 samples per scan
    scan_rate = 8  # single scan, seconds
    scan_angle = -48.3  # swath, degrees
    sampling_interval = 0.2  # single view, seconds
    sync_time = 0.00355  # delay before the actual scan starts

    if scan_points is None:
        scan_points = np.arange(0, scan_len)

    # build the instrument (scan angles)
    samples = np.vstack(((scan_points / (scan_len * 0.5 - 0.5) - 1)
                         * np.deg2rad(scan_angle),
                         np.zeros((len(scan_points),))))
    samples = np.tile(samples[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

    # building the corresponding times array
    offset = np.arange(scans_nb) * scan_rate
    times = (np.tile(scan_points * sampling_interval + sync_time,
                     [np.int32(scans_nb), 1])
             + np.expand_dims(offset, 1))

    # build the scan geometry object
    return ScanGeometry(samples, times)


################################################################
#
#   MHS
#
################################################################

def mhs(scans_nb, scan_points=None):
    """ Describe MHS instrument geometry

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

    scan_len = 90  # 90 samples per scan
    scan_rate = 8 / 3.  # single scan, seconds
    scan_angle = -49.444  # swath, degrees
    sampling_interval = (8 / 3. - 1) / 90.  # single view, seconds
    sync_time = 0.0  # delay before the actual scan starts - don't know! FIXME!

    if scan_points is None:
        scan_points = np.arange(0, scan_len)

    # build the instrument (scan angles)
    samples = np.vstack(((scan_points / (scan_len * 0.5 - 0.5) - 1) * np.deg2rad(scan_angle),
                         np.zeros((len(scan_points),))))
    samples = np.tile(samples[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

    # building the corresponding times array
    offset = np.arange(scans_nb) * scan_rate
    times = (np.tile(scan_points * sampling_interval + sync_time, [np.int32(scans_nb), 1]) + np.expand_dims(offset, 1))

    # scan_angles = np.linspace(-np.deg2rad(scan_angle), np.deg2rad(scan_angle), scan_len)[scan_points]

    # samples = np.vstack((scan_angles, np.zeros(len(scan_points) * 1,)))
    # samples = np.tile(samples[:, np.newaxis, :], [1, np.int(scans_nb), 1])

    # # building the corresponding times array
    # offset = np.arange(scans_nb) * scan_rate
    # times = (np.tile(scan_points * sampling_interval, [np.int(scans_nb), 1])
    #          + np.expand_dims(offset, 1))

    # build the scan geometry object
    return ScanGeometry(samples, times)


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

    scan_len = 56  # 56 samples per scan
    scan_rate = 6.4  # single scan, seconds
    scan_angle = -49.5  # swath, degrees
    sampling_interval = abs(scan_rate) / scan_len  # single view, seconds

    if scan_points is None:
        scan_points = np.arange(0, scan_len)

    # build the instrument (scan angles)
    samples = np.vstack(((scan_points / (scan_len * 0.5 - 0.5) - 1)
                         * np.deg2rad(scan_angle),
                         np.zeros((len(scan_points),))))
    samples = np.tile(samples[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

    # building the corresponding times array
    offset = np.arange(scans_nb) * scan_rate
    times = (np.tile(scan_points * sampling_interval, [np.int32(scans_nb), 1])
             + np.expand_dims(offset, 1))

    # build the scan geometry object
    return ScanGeometry(samples, times)


################################################################
#
#   ATMS
#
################################################################

def atms(scans_nb, scan_points=None):
    """ Describe ATMS instrument geometry
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

    scan_len = 96  # 96 samples per scan
    scan_rate = 8 / 3.  # single scan, seconds
    scan_angle = -52.7  # swath, degrees
    sampling_interval = 18e-3  # single view, seconds

    if scan_points is None:
        scan_points = np.arange(0, scan_len)

    # build the instrument (scan angles)
    scan_angles = np.linspace(-np.deg2rad(scan_angle), np.deg2rad(scan_angle), scan_len)[scan_points]

    samples = np.vstack((scan_angles, np.zeros(len(scan_points) * 1,)))
    samples = np.tile(samples[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

    # building the corresponding times array
    offset = np.arange(scans_nb) * scan_rate
    times = (np.tile(scan_points * sampling_interval, [np.int32(scans_nb), 1])
             + np.expand_dims(offset, 1))

    # build the scan geometry object
    return ScanGeometry(samples, times)


################################################################
#
#   MWHS-2
#
################################################################

def mwhs2(scans_nb, scan_points=None):
    """Describe MWHS-2 instrument geometry

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

    scan_len = 98  # 98 samples per scan
    scan_rate = 8 / 3.  # single scan, seconds
    scan_angle = -53.35  # swath, degrees
    sampling_interval = (8 / 3. - 1) / 98.  # single view, seconds
    # sampling_interval = 17.449e-3  # single view, seconds
    sync_time = 0.0  # delay before the actual scan starts - don't know! FIXME!

    if scan_points is None:
        scan_points = np.arange(0, scan_len)

    # build the instrument (scan angles)
    samples = np.vstack(((scan_points / (scan_len * 0.5 - 0.5) - 1)
                         * np.deg2rad(scan_angle),
                         np.zeros((len(scan_points),))))
    samples = np.tile(samples[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

    # building the corresponding times array
    offset = np.arange(scans_nb) * scan_rate
    times = (np.tile(scan_points * sampling_interval + sync_time,
                     [np.int32(scans_nb), 1])
             + np.expand_dims(offset, 1))

    # # build the instrument (scan angles)
    # scan_angles = np.linspace(-np.deg2rad(scan_angle), np.deg2rad(scan_angle), scan_len)[scan_points]

    # samples = np.vstack((scan_angles, np.zeros(len(scan_points) * 1,)))
    # samples = np.tile(samples[:, np.newaxis, :], [1, np.int(scans_nb), 1])

    # # building the corresponding times array
    # offset = np.arange(scans_nb) * scan_rate
    # times = (np.tile(scan_points * sampling_interval, [np.int(scans_nb), 1])
    #          + np.expand_dims(offset, 1))

    # build the scan geometry object
    return ScanGeometry(samples, times)


################################################################
#
#   OLCI
#
################################################################


def olci(scans_nb, scan_points=None):
    """Definition of the OLCI instrument.

    Source: Sentinel-3 OLCI Coverage
    https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-olci/coverage
    """

    if scan_points is None:
        scan_len = 4000  # samples per scan
        scan_points = np.arange(4000)
    else:
        scan_len = len(scan_points)
    # scan_rate = 0.044  # single scan, seconds
    scan_angle_west = 46.5  # swath, degrees
    scan_angle_east = -22.1  # swath, degrees
    # sampling_interval = 18e-3  # single view, seconds
    # build the olci instrument scan line angles
    scanline_angles = np.linspace(np.deg2rad(scan_angle_west),
                                  np.deg2rad(scan_angle_east), scan_len)
    inst = np.vstack((scanline_angles, np.zeros(scan_len,)))

    inst = np.tile(inst[:, np.newaxis, :], [1, np.int32(scans_nb), 1])

    # building the corresponding times array
    # times = (np.tile(scan_points * 0.000025 + 0.0025415, [scans_nb, 1])
    #         + np.expand_dims(offset, 1))

    times = np.tile(np.zeros_like(scanline_angles), [np.int32(scans_nb), 1])
    # if apply_offset:
    #     offset = np.arange(np.int(scans_nb)) * frequency
    #     times += np.expand_dims(offset, 1)

    return ScanGeometry(inst, times)


def ascat(scan_nb, scan_points=None):
    """ASCAT make two scans one to the left and one to the right of the
    sub-satellite track.

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
