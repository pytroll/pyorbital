#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011 SMHI

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

"""Astronomy module.
Parts taken from http://www.geoastro.de/elevaz/basics/index.htm
"""

import datetime
import numpy as np

def jdays2000(current_time):
    """Get the days since 2000.
    """
    return _days(current_time - datetime.datetime(2000, 1, 1, 12, 0))
    

def jdays(current_time):
    """Get the julian day of *current_time*.
    """
    return jdays2000(current_time) + 2451545

def _days(dt):
    """Get the days (floating point) from *d_t*.
    """
    return (dt.days +
            (dt.seconds +
             dt.microseconds / (1000000.0)) / (24 * 3600.0))

def gmst(current_time):
    """Greenwich mean sidereal current_time, in radians.
    http://celestrak.com/columns/v02n02/
    """
    now = current_time
    #now = datetime.datetime(1995, 10, 1, 9, 0)
    now0 = datetime.datetime(now.year, now.month, now.day)
    epoch = datetime.datetime(2000, 1, 1, 12, 0)
    du2 = _days(now - epoch)
    d_u = _days(now0 - epoch)

    dus = (du2 - d_u) * 86400
    t_u = d_u / 36525.0
    theta_g_0 = (24110.54841 + t_u * (8640184.812866 +
                                      t_u * (0.093104 - t_u * 6.2 * 10e-6)))
    theta_g = (theta_g_0 + dus * 1.00273790934) % 86400
    return (theta_g / 86400.0) * 2 * np.pi

def lmst(current_time, longitude):
    """Local mean sidereal time, computed from *current_time* and *longitude*.
    In radians.
    """
    return gmst(current_time) + longitude


def sun_ecliptic_longitude(current_time):
    """Ecliptic longitude of the sun at *current_time*.
    """
    jdate = jdays2000(current_time) / 36525.0
    # mean anomaly, rad
    m_a = np.deg2rad(357.52910 +
                     35999.05030*jdate -
                     0.0001559*jdate*jdate -
                     0.00000048*jdate*jdate*jdate)
    # mean longitude, deg
    l_0 = 280.46645 + 36000.76983*jdate + 0.0003032*jdate*jdate
    d_l = ((1.914600 - 0.004817*jdate - 0.000014*jdate*jdate)*np.sin(m_a) +
           (0.019993 - 0.000101*jdate)*np.sin(2*m_a) + 0.000290*np.sin(3*m_a))
    # true longitude, deg
    l__ = l_0 + d_l
    return np.deg2rad(l__)

def sun_ra_dec(current_time):
    """Right ascension and declination of the sun at *current_time*.
    """
    jdate = jdays2000(current_time) / 36525.0
    eps = np.deg2rad(23.0 + 26.0/60.0 + 21.448/3600.0 -
                     (46.8150*jdate + 0.00059*jdate*jdate -
                      0.001813*jdate*jdate*jdate) / 3600)
    eclon = sun_ecliptic_longitude(current_time)
    x__ = np.cos(eclon)
    y__ = np.cos(eps) * np.sin(eclon)
    z__ = np.sin(eps) * np.sin(eclon)
    r__ = np.sqrt(1.0 - z__ * z__)
    # sun declination
    declination = np.arctan2(z__, r__)
    # right ascension
    right_ascension = 2 * np.arctan2(y__, (x__ + r__))
    return right_ascension, declination

def local_hour_angle(current_time, longitude, right_ascension):
    """Hour angle at *current_time* for the given *longitude* and
    *right_ascension*
    """
    return lmst(current_time, longitude) - right_ascension

def get_alt_az(current_time, lon, lat):
    """Return sun altitude and azimuth from *current_time*, *lon*, and *lat*.
    """
    ra_, dec = sun_ra_dec(current_time)
    h__ = local_hour_angle(current_time, lon, ra_)
    return (np.arcsin(np.sin(lat)*np.sin(dec) +
                      np.cos(lat) * np.cos(dec) * np.cos(h__)),
            np.arctan2(-np.sin(h__), (np.cos(lat)*np.tan(dec) -
                                    np.sin(lat)*np.cos(h__))))

def cos_zen(current_time, lon, lat):
    """Cosine of the sun-zenith angle for *lon*, *lat* at *current_time*.
    """
    r_a, dec = sun_ra_dec(current_time)
    h__ = local_hour_angle(current_time, lon, r_a)
    return (np.sin(lat)*np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h__))

def observer_position(time, lon, lat, alt):
    """Calculate observer ECI position.
    http://celestrak.com/columns/v02n02/
    """
    theta = (gmst(time) + lon)%(2*np.pi)
    c = 1/np.sqrt(1 + F*(F-2)*np.sin(lat)**2)
    sq = c*(1 - F)**2

    achcp = (XKMPER*c + alt)*np.cos(lat)
    x = achcp*np.cos(theta)  # kilometers
    y = achcp*np.sin(theta)
    z = (XKMPER*sq + alt)*np.sin(lat)

    vx = -MFACTOR*y  # kilometers/second
    vy = MFACTOR*x
    vz = 0

    return (x, y, z), (vx, vy, vz)
    
        
if __name__ == '__main__':
    pass
