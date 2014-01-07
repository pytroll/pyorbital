#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011, 2012, 2013, 2014.

# Author(s):

#   Esben S. Nielsen <esn@dmi.dk>
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


import logging
import datetime
import urllib2

tle_urls = ('http://celestrak.com/NORAD/elements/weather.txt',
            'http://celestrak.com/NORAD/elements/resource.txt')

logger = logging.getLogger(__name__)

def read(platform, tle_file=None, line1=None, line2=None):
    """Read TLE for *satellite* from *tle_file*, from *line1* and *line2*, or
    from internet if none is provided.
    """
    return Tle(platform, tle_file=tle_file, line1=line1, line2=line2)

def fetch(destination):
    """fetch TLE from internet and save it to *destination*.
    """
    with open(destination, "w") as dest:
        for url in tle_urls:
            response = urllib2.urlopen(url)
            dest.write(response.read())

class ChecksumError(Exception):
    pass


class Tle(object):
    """Class holding TLE objects.
    """    

    def __init__(self, platform, tle_file=None, line1=None, line2=None):
        platform = platform.strip().upper()

        if line1 is not None and line2 is not None:
            tle = line1.strip() + "\n" + line2.strip()
        else:
            if tle_file:
                urls = (tle_file,)
                open_func = open
            else:
                logger.debug("Fetch tle from the internet.")
                urls = tle_urls
                open_func = urllib2.urlopen
            
            tle = ""
            for url in urls:
                fp = open_func(url)
                for l0 in fp:
                    l1, l2 = fp.next(), fp.next()
                    if l0.strip() == platform:
                        tle = l1.strip() + "\n" + l2.strip()
                        break
                fp.close()
                if tle:
                    break
            
            if not tle:
                raise AttributeError, "Found no TLE entry for '%s'" % platform

        self._platform = platform
        self._line1, self._line2 = tle.split('\n')
        self._checksum()
        self._read_tle()

    @property
    def line1(self):
        return self._line1

    @property
    def line2(self):
        return self._line2

    @property
    def platform(self):
        return self._platform

    def _checksum(self):
        """Performs the checksum for the current TLE.
        """
        for line in [self.line1, self.line2]:
            check = 0
            for char in line[:-1]:
                if char.isdigit():
                    check += int(char)
                if char == "-":
                    check += 1

            if (check % 10) != int(line[-1]):
                raise ChecksumError(self._platform + " " + line)

    def _read_tle(self):

        def _read_tle_decimal(rep):
            if rep[0] in ["-", " "]:
                val = rep[0] + "." + rep[1:-2] + "e" + rep[-2:]
            else:
                val = "." + rep[:-2] + "e" + rep[-2:]

            return float(val)
        
        self.satnumber = self._line1[2:7]
        self.classification = self._line1[7]
        self.id_launch_year = self._line1[9:11]
        self.id_launch_number = self._line1[11:14]
        self.id_launch_piece = self._line1[14:17]
        self.epoch_year = self._line1[18:20]
        self.epoch_day = float(self._line1[20:32])
        self.epoch = (datetime.datetime.strptime(self.epoch_year, "%y") +
                             datetime.timedelta(days=self.epoch_day - 1))
        self.mean_motion_derivative = float(self._line1[33:43])
        self.mean_motion_sec_derivative = _read_tle_decimal(self._line1[44:52])
        self.bstar = float(self._line1[53] + "." + self._line1[54:59] + "e" + self._line1[59:61])
        _read_tle_decimal(self._line1[53:61])
        try:
            self.ephemeris_type = int(self._line1[62])
        except ValueError:
            self.ephemeris_type = 0
        self.element_number = int(self._line1[64:68])

        self.inclination = float(self._line2[8:16])
        self.right_ascension = float(self._line2[17:25])
        self.excentricity = int(self._line2[26:33]) * 10 ** -7
        self.arg_perigee = float(self._line2[34:42])
        self.mean_anomaly = float(self._line2[43:51])
        self.mean_motion = float(self._line2[52:63])
        self.orbit = int(self._line2[63:68])

    def __str__(self):
        import pprint, StringIO
        s = StringIO.StringIO()
        d = dict(([(k, v) for k, v in self.__dict__.items() if k[0] != '_']))
        pprint.pprint(d, s)
        return s.getvalue()[:-1]

if __name__ == '__main__':
    tle = read('noaa 19')
    print tle
