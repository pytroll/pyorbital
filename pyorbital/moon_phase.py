#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2016 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>

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

"""Calculations of the phase of the moon
"""

from numpy import deg2rad, sin, cos, tan, arctan, fabs
import numpy as np
import datetime
import collections

import math
PI = math.pi
SMALL_FLOAT = 1.0e-12


def seqJulian(dobj):
    """Returns the number of julian days for the specified date-time.
    Input = datetime object.
    """

    if isinstance(dobj, datetime.datetime):
        dobj = np.array([dobj], 'datetime64[us]')
    elif isinstance(dobj, np.ndarray):
        dobj = dobj.astype('datetime64[us]')
    elif isinstance(dobj, collections.Sequence):
        dobj = np.array(dobj, 'datetime64[us]')

    day = (
        dobj - dobj.astype('datetime64[M]')).astype('f') / (24. * 3600. * 1000) + 1.0
    month = (
        dobj.astype('datetime64[M]') - dobj.astype('datetime64[Y]')).astype('f') + 1
    year = dobj.astype('datetime64[Y]').astype('f') + 1970

    year = np.where(np.less(month, 3), year - 1, year)
    month = np.where(np.less(month, 3), month + 12, month)

    cond1 = np.less(year, 1582)
    cond2 = np.logical_and(np.equal(year, 1582), np.greater(month, 10))
    cond3 = np.logical_and(np.logical_and(np.equal(year, 1582), np.equal(month, 10)),
                           np.greater(day, 15))

    a__ = np.divide(year, 100).astype('i')
    b__ = np.where(np.logical_or(np.logical_or(cond1, cond2), cond3),
                   2 - a__ + a__ / 4, -10)

    c__ = (365.25 * year).astype('i')
    e__ = (30.6001 * (month + 1)).astype('i')
    return b__ + c__ + e__ + day + 1720994.5


def Julian(dobj):
    """Returns the number of julian days for the specified date-time.
    Input = datetime object.
    """

    year = dobj.year
    month = dobj.month
    day = dobj.day + (dobj.hour / 24. + dobj.minute / (24 * 60.) +
                      dobj.second / (24 * 3600.) +
                      dobj.microsecond / (24 * 3600 * 1000000.))

    if month < 3:
        year = year - 1
        month += 12

    if (year > 1582 or
            (year == 1582 and month > 10) or
            (year == 1582 and month == 10 and day > 15)):
        a__ = int(year / 100)
        b__ = 2 - a__ + a__ / 4

    c__ = int(365.25 * year)
    e__ = int(30.6001 * (month + 1))
    return b__ + c__ + e__ + day + 1720994.5


def sun_position(jday):
    """Get sun position"""

    # double n,x,e,l,dl,v;
    # double m2;
    # int i;

    if isinstance(jday, collections.Sequence):
        jday = np.array(jday)
    elif not isinstance(jday, np.ndarray):
        jday = np.array([jday], 'f')

    n__ = 360. / 365.2422 * jday
    i__ = np.divide(n__, 360.).astype('i')
    n__ = n__ - i__ * 360.0
    x__ = n__ - 3.762863
    x__ = np.where(np.less(x__, 0), x__ + 360., x__)

    x__ = deg2rad(x__)
    e__ = x__
    while 1:
        dl_ = e__ - .016718 * sin(e__) - x__
        e__ = e__ - dl_ / (1 - .016718 * cos(e__))
        if np.alltrue(np.less(np.fabs(dl_), SMALL_FLOAT)):
            break

    v__ = 360. / PI * arctan(1.01686011182 * tan(e__ / 2))
    sunpos = v__ + 282.596403
    i__ = np.divide(sunpos, 360.).astype('i')
    sunpos = sunpos - i__ * 360.0

    if isinstance(sunpos, np.ndarray) and len(sunpos) == 1:
        return sunpos[0]
    else:
        return sunpos

    return sunpos

# def moon_position(jday, lsun):
#     """Get the moon position"""

#     ms_ = 0.985647332099 * jday - 3.762863
#     if ms_ < 0:
#         ms_ = ms_ + 360.0

#     mpos = 13.176396 * jday + 64.975464
#     i__ = int(mpos / 360.0)
#     mpos = mpos - i__ * 360.0
#     if mpos < 0:
#         mpos = mpos + 360.0

#     mm_ = mpos - 0.1114041 * jday - 349.383063
#     i__ = int(mm_ / 360.0)
#     mm_ = mm_ - i__ * 360.0

#     ev_ = 1.2739 * sin(deg2rad(2 * (mpos - lsun) - mm_))
#     sms = sin(deg2rad(ms_))
#     ae_ = 0.1858 * sms
#     mm_ = mm_ + ev_ - ae_ - 0.37 * sms
#     ec_ = 6.2886 * sin(deg2rad(mm_))
#     mpos = mpos + ev_ + ec_ - ae_ + 0.214 * sin(deg2rad(2 * mm_))
#     mpos = 0.6583 * sin(deg2rad(2 * (mpos - lsun))) + mpos

#     return mpos


def moon_position(jday, lsun):
    """Get the moon position"""

    if isinstance(jday, collections.Sequence):
        jday = np.array(jday)
    elif not isinstance(jday, np.ndarray):
        jday = np.array([jday], 'f')

    if isinstance(lsun, collections.Sequence):
        lsun = np.array(lsun)
    elif not isinstance(lsun, np.ndarray):
        lsun = np.array([lsun], 'f')

    ms_ = 0.985647332099 * jday - 3.762863
    ms_ = np.where(np.less(ms_, 0), ms_ + 360.0, ms_)

    mpos = 13.176396 * jday + 64.975464
    i__ = np.divide(mpos, 360.0).astype('i')
    mpos = mpos - i__ * 360.0
    mpos = np.where(np.less(mpos, 0), mpos + 360.0, mpos)

    mm_ = mpos - 0.1114041 * jday - 349.383063
    i__ = np.divide(mm_, 360.0).astype('i')
    mm_ = mm_ - i__ * 360.0

    ev_ = 1.2739 * sin(deg2rad(2 * (mpos - lsun) - mm_))
    sms = sin(deg2rad(ms_))
    ae_ = 0.1858 * sms
    mm_ = mm_ + ev_ - ae_ - 0.37 * sms
    ec_ = 6.2886 * sin(deg2rad(mm_))
    mpos = mpos + ev_ + ec_ - ae_ + 0.214 * sin(deg2rad(2 * mm_))
    mpos = 0.6583 * sin(deg2rad(2 * (mpos - lsun))) + mpos

    if isinstance(mpos, np.ndarray) and len(mpos) == 1:
        return mpos[0]
    else:
        return mpos


def moon_phase(dobj):
    """Calculate the phase of the moon"""

    jday = Julian(dobj) - 2444238.5
    lsun = sun_position(jday)
    lmoon = moon_position(jday, lsun)

    phase = (1.0 - cos(deg2rad(lmoon - lsun))) / 2.0
    if isinstance(phase, np.ndarray) and len(phase) == 1:
        phase = phase[0]
    return phase
