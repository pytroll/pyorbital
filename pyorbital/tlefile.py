#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011, 2012, 2013, 2014, 2015.

# Author(s):

#   Esben S. Nielsen <esn@dmi.dk>
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Panu Lahtinen <panu.lahtinen@fmi.fi>

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


import io
import logging
import datetime
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
import os
import glob

TLE_URLS = ('http://celestrak.com/NORAD/elements/weather.txt',
            'http://celestrak.com/NORAD/elements/resource.txt')

LOGGER = logging.getLogger(__name__)


def read_platform_numbers(in_upper=False, num_as_int=False):
    '''Read platform numbers from $PPP_CONFIG_DIR/platforms.txt if available.
    '''

    out_dict = {}
    if "PPP_CONFIG_DIR" in os.environ:
        platform_file = os.path.join(os.environ["PPP_CONFIG_DIR"],
                                     "platforms.txt")
        try:
            fid = open(platform_file, 'r')
        except IOError:
            LOGGER.error("Platform file %s not found.", platform_file)
            return out_dict
        for row in fid:
            # skip comment lines
            if not row.startswith('#'):
                parts = row.split()
                if in_upper:
                    parts[0] = parts[0].upper()
                if num_as_int:
                    parts[1] = int(parts[1])
                out_dict[parts[0]] = parts[1]
        fid.close()

    return out_dict


SATELLITES = read_platform_numbers(in_upper=True, num_as_int=False)
'''
The platform numbers are given in a file $PPP_CONFIG/platforms.txt
in the following format:

# Mappings between satellite catalogue numbers and corresponding
# platform names from OSCAR.
ALOS-2 39766
CloudSat 29107
CryoSat-2 36508
CSK-1 31598
CSK-2 32376
CSK-3 33412
CSK-4 37216
DMSP-F15 25991
DMSP-F16 28054
DMSP-F17 29522
DMSP-F18 35951
DMSP-F19 39630
EOS-Aqua 27424
EOS-Aura 28376
EOS-Terra 25994
FY-2D 29640
FY-2E 33463
FY-2F 38049
FY-2G 40367
FY-3A 32958
FY-3B 37214
FY-3C 39260
GOES-13 29155
GOES-14 35491
GOES-15 36411
Himawari-6 28622
Himawari-7 28937
Himawari-8 40267
INSAT-3A 27714
INSAT-3C 27298
INSAT-3D 39216
JASON-2 33105
Kalpana-1 27525
Landsat-7 25682
Landsat-8 39084
Meteosat-7 24932
Meteosat-8 27509
Meteosat-9 28912
Meteosat-10 38552
Metop-A 29499
Metop-B 38771
NOAA-15 25338
NOAA-16 26536
NOAA-17 27453
NOAA-18 28654
NOAA-19 33591
RadarSat-2 32382
Sentinel-1A 39634
SMOS 36036
SPOT-5 27421
SPOT-6 38755
SPOT-7 40053
Suomi-NPP 37849
TanDEM-X 36605
TerraSAR-X 31698
'''


def read(platform, tle_file=None, line1=None, line2=None):
    """Read TLE for *satellite* from *tle_file*, from *line1* and *line2*, from
   the newest file provided in the TLES pattern, or from internet if none is
   provided.
   """
    return Tle(platform, tle_file=tle_file, line1=line1, line2=line2)


def fetch(destination):
    """fetch TLE from internet and save it to *destination*.
   """
    with open(destination, "w") as dest:
        for url in TLE_URLS:
            response = urlopen(url)
            dest.write(response.read())


class ChecksumError(Exception):
    '''ChecksumError.
    '''
    pass


class Tle(object):
    """Class holding TLE objects.
   """

    def __init__(self, platform, tle_file=None, line1=None, line2=None):
        self._platform = platform.strip().upper()
        self._tle_file = tle_file
        self._line1 = line1
        self._line2 = line2

        self.satnumber = None
        self.classification = None
        self.id_launch_year = None
        self.id_launch_number = None
        self.id_launch_piece = None
        self.epoch_year = None
        self.epoch_day = None
        self.epoch = None
        self.mean_motion_derivative = None
        self.mean_motion_sec_derivative = None
        self.bstar = None
        self.ephemeris_type = None
        self.element_number = None
        self.inclination = None
        self.right_ascension = None
        self.excentricity = None
        self.arg_perigee = None
        self.mean_anomaly = None
        self.mean_motion = None
        self.orbit = None

        self._read_tle()
        self._checksum()
        self._parse_tle()

    @property
    def line1(self):
        '''Return first TLE line.'''
        return self._line1

    @property
    def line2(self):
        '''Return second TLE line.'''
        return self._line2

    @property
    def platform(self):
        '''Return satellite platform name.'''
        return self._platform

    def _checksum(self):
        """Performs the checksum for the current TLE.
        """
        for line in [self._line1, self._line2]:
            check = 0
            for char in line[:-1]:
                if char.isdigit():
                    check += int(char)
                if char == "-":
                    check += 1

            if (check % 10) != int(line[-1]):
                raise ChecksumError(self._platform + " " + line)

    def _read_tle(self):
        '''Read TLE data.
        '''

        if self._line1 is not None and self._line2 is not None:
            tle = self._line1.strip() + "\n" + self._line2.strip()
        else:
            def _open(filename):
                return io.open(filename, 'rb')

            if self._tle_file:
                urls = (self._tle_file,)
                open_func = _open
            elif "TLES" in os.environ:
                # TODO: get the TLE file closest in time to the actual satellite
                # overpass, NOT the latest!
                urls = (max(glob.glob(os.environ["TLES"]),
                            key=os.path.getctime), )
                LOGGER.debug("Reading TLE from %s", urls[0])
                open_func = _open
            else:
                LOGGER.debug("Fetch TLE from the internet.")
                urls = TLE_URLS
                open_func = urlopen

            tle = ""
            designator = "1 " + SATELLITES.get(self._platform, '')
            for url in urls:
                fid = open_func(url)
                for l_0 in fid:
                    l_0 = l_0.decode('utf-8')
                    if l_0.strip() == self._platform:
                        l_1 = next(fid).decode('utf-8')
                        l_2 = next(fid).decode('utf-8')
                        tle = l_1.strip() + "\n" + l_2.strip()
                        break
                    if(self._platform in SATELLITES and
                       l_0.strip().startswith(designator)):
                        l_1 = l_0
                        l_2 = next(fid).decode('utf-8')
                        tle = l_1.strip() + "\n" + l_2.strip()
                        LOGGER.debug("Found platform %s, ID: %s",
                                     self._platform,
                                     SATELLITES[self._platform])
                        break
                fid.close()
                if tle:
                    break

            if not tle:
                raise KeyError("Found no TLE entry for '%s'" % self._platform)

        self._line1, self._line2 = tle.split('\n')

    def _parse_tle(self):
        '''Parse values from TLE data.
        '''
        def _read_tle_decimal(rep):
            '''Convert *rep* to decimal value.
            '''
            if rep[0] in ["-", " ", "+"]:
                digits = rep[1:-2].strip()
                val = rep[0] + "." + digits + "e" + rep[-2:]
            else:
                digits = rep[:-2].strip()
                val = "." + digits + "e" + rep[-2:]

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
        self.bstar = _read_tle_decimal(self._line1[53:61])
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
        import pprint
        import sys
        if sys.version_info < (3, 0):
            from StringIO import StringIO
        else:
            from io import StringIO
        s_var = StringIO()
        d_var = dict(([(k, v) for k, v in
                       list(self.__dict__.items()) if k[0] != '_']))
        pprint.pprint(d_var, s_var)
        return s_var.getvalue()[:-1]

def main():
    '''Main for testing TLE reading.
    '''
    tle_data = read('Noaa-19')
    print(tle_data)

if __name__ == '__main__':
    main()
