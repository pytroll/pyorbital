#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Mikhail Itkin <itkin.m@gmail.com>

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

# number of instrument scans to use.
scans_nb = 10

################################################################
#
#   AVHRR
#
################################################################

def avhrr(scans_nb, scan_points, scan_angle=55.37, decimate=1):
    """Definition of the avhrr instrument.

    Source: NOAA KLM User's Guide, Appendix J
    http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/j/app-j.htm
    """
    # build the avhrr instrument (scan angles)
    avhrr_inst = np.vstack(((scan_points / 1023.5 - 1)
                            * np.deg2rad(-scan_angle),
                            np.zeros((len(scan_points),)))).transpose()
    avhrr_inst = np.tile(avhrr_inst, [scans_nb, 1])

    # building the corresponding times array
    offset = np.arange(scans_nb) * decimate / 6.0
    #times = (np.tile(scan_points * 0.000025 + 0.0025415, [scans_nb, 1])
    #         + np.expand_dims(offset, 1))
    times = (np.tile(scan_points * 0.000025, [scans_nb, 1])
             + np.expand_dims(offset, 1))

    return ScanGeometry(avhrr_inst, times.ravel())

################################################################
#### avhrr, all pixels

# we take all pixels
scan_points = np.arange(2048)

# build the scan geometry object
avhrr_all_geom = avhrr(scans_nb, scan_points)

################################################################
#### avhrr, edge pixels

# we take only edge pixels
scan_points = np.array([0, 2047])


# build the scan geometry object
avhrr_edge_geom = avhrr(scans_nb, scan_points)

################################################################
#### avhrr, every 40th pixel from the 24th (aapp style)

# we take only every 40th pixel
scan_points = np.arange(24, 2048, 40)

# build the scan geometry object
avhrr_40_geom = avhrr(scans_nb, scan_points)

################################################################
#
#   VIIRS
#
################################################################

def viirs(scans_nb, scan_indices=slice(0, None)):
    """Describe VIIRS instrument geometry, I-band.

    """

    entire_width = np.arange(6400)
    scan_points = entire_width[scan_indices]
    
    across_track = (scan_points / 3199.5 - 1) * np.deg2rad(-55.84)
    y_max_angle = np.arctan2(11.87/2, 824.0)
    along_track = np.array([-y_max_angle, 0, y_max_angle])

    scan_pixels = len(scan_points)

    scan = np.vstack((np.tile(across_track, scan_pixels),
                      np.repeat(along_track, 6400))).T
    
    npp = np.tile(scan, [scans_nb, 1])
    
    # from the timestamp in the filenames, a granule takes 1:25.400 to record
    # (85.4 seconds) so 1.779166667 would be the duration of 1 scanline
    # dividing the duration of a single scan by a width of 6400 pixels results
    # in 0.0002779947917 seconds for each column of 32 pixels in the scanline

    # the individual times per pixel are probably wrong, unless the scanning
    # behaves the same as for AVHRR, The VIIRS sensor rotates to allow internal
    # calibration before each scanline. This would imply that the scanline
    # always moves in the same direction.  more info @
    # http://www.eoportal.org/directory/pres_NPOESSNationalPolarorbitingOperationalEnvironmentalSatelliteSystem.html

    offset = np.arange(scans_nb) * 1.779166667
    times = (np.tile(scan_points * 0.0002779947917, [scans_nb, scan_pixels])
             + np.expand_dims(offset, 1))

    # build the scan geometry object
    return ScanGeometry(npp, times.ravel())


################################################################
#
#   AMSU-A
#
################################################################

def amsua(scans_nb, edges_only=False):
    """ Describe AMSU-A instrument geometry
    
    Parameters:
       scans_nb | int -  number of scan lines
     
     Keywords:
     * edges_only - use only edge pixels

    Returns:
       pyorbital.geoloc.ScanGeometry object
    
    """

    scan_len  = 30 # 30 samples per scan
    scan_rate = 8 # single scan, seconds
    scan_angle = -48.3 # swath, degrees
    sampling_interval = 0.2 # single view, seconds
    sync_time = 0.00355 # delay before the actual scan starts

    if edges_only:
        scan_points = np.array([0, scan_len - 1])
    else:
        scan_points = np.arange(0, scan_len)

    # build the instrument (scan angles)
    samples = np.vstack(((scan_points / (scan_len*0.5-0.5) - 1)
                         * np.deg2rad(scan_angle),
                         np.zeros((len(scan_points),)))).transpose()
    samples = np.tile(samples, [scans_nb, 1])

    # building the corresponding times array
    offset = np.arange(scans_nb) * scan_rate
    times = (np.tile(scan_points * sampling_interval + sync_time, [scans_nb, 1])
	         + np.expand_dims(offset, 1))

    # build the scan geometry object
    return ScanGeometry(samples, times.ravel())
