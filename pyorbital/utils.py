#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2016 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c20671.ad.smhi.se>

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

"""
"""


def JulianToDate(jday):
    """Get the current datetime object from the julian day"""
    import datetime

    jday = jday + 0.5
    jdi = long(jday)
    if (jdi > 2299160):
        a__ = long((jdi - 1867216.25) / 36524.25)
        b__ = long(jdi + 1 + a__ - a__ / 4)
    else:
        b__ = jdi

    c__ = long(b__ + 1524)
    d__ = long((c__ - 122.1) / 365.25)
    e__ = long(365.25 * d__)
    g__ = long((c__ - e__) / 30.6001)
    g1_ = long(30.6001 * g__)

    day = c__ - e__ - g1_
    fhour = float((jday - jdi) * 24.0)
    if (g__ <= 13):
        month = g__ - 1
    else:
        month = g__ - 13

    if (month > 2):
        year = d__ - 4716
    else:
        year = d__ - 4715

    hour = int(fhour)
    minutes_of_hour = 60 * (fhour - int(hour))
    minutes = int(minutes_of_hour)
    seconds = 60 * (minutes_of_hour - minutes)
    microsec = int((seconds - int(seconds)) * 1000000)
    seconds = int(seconds)

    return datetime.datetime(year, month, day,
                             hour=hour, minute=minutes,
                             second=seconds, microsecond=microsec)
