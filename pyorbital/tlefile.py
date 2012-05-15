#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011, 2012.

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


import datetime

tle_urls = ('http://celestrak.com/NORAD/elements/weather.txt',
            'http://celestrak.com/NORAD/elements/resource.txt')

def read(platform, tle_file=None, line1=None, line2=None):
    """Read TLE for *satellite* from *tle_file*, from *line1* and *line2*, or
    from internet if none is provided.
    """
    return Tle(platform, tle_file=tle_file, line1=line1, line2=line2)

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
                import urllib2
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

    def _read_tle(self):

        def _read_tle_decimal(rep):
            num = int(rep[:-2]) * 1.0e-5
            exp = int(rep[-2:])
            return num * 10 ** exp

        lines = self._line1.split() + self._line2.split()
        if len(lines) < 17:
            raise IOError("To few items in the TLE entry, found '%d', expected at least '17'"%len(lines)) 

        self.satnumber = lines[1][:5]
        self.classification = lines[1][5:]
        self.id_launch_year = lines[2][:2]
        self.id_launch_number = lines[2][2:5]
        self.id_launch_piece = lines[2][5:]
        self.epoch_year = int(lines[3][:2])
        self.epoch_day = float(lines[3][2:])
        self.epoch = (datetime.datetime.strptime(lines[3][:2], "%y") +
                             datetime.timedelta(days=float(lines[3][2:]) - 1))
        self.mean_motion_derivative = float(lines[4])
        self.mean_motion_sec_derivative = _read_tle_decimal(lines[5])
        self.bstar = _read_tle_decimal(lines[6])
        self.ephemeris_type = int(lines[7])
        self.element_number = int(lines[8][:-1])

        self.inclination = float(lines[11])
        self.right_ascension = float(lines[12])
        self.excentricity = int(lines[13]) * 10 ** -7
        self.arg_perigee = float(lines[14])
        self.mean_anomaly = float(lines[15])
        if len(lines) > 17:
            self.orbit = int(lines[17][:-1])
            self.mean_motion = float(lines[16])
        else:
            self.orbit = int(lines[16][11:-1])
            self.mean_motion = float(lines[16][:11])

    def __str__(self):
        import pprint, StringIO
        s = StringIO.StringIO()
        d = dict(([(k, v) for k, v in self.__dict__.items() if k[0] != '_']))
        pprint.pprint(d, s)
        return s.getvalue()[:-1]

if __name__ == '__main__':
    tle = read('noaa 19')
    print tle
