#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2011 - 2018
#
# Author(s):
#
#   Esben S. Nielsen <esn@dmi.dk>
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Panu Lahtinen <panu.lahtinen@fmi.fi>
#   Will Evonosky <william.evonosky@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Classes and functions for handling TLE files."""

import io
import logging
import datetime as dt
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
import os
import glob
import numpy as np
import requests
import sqlite3

TLE_URLS = ('http://www.celestrak.com/NORAD/elements/active.txt',
            'http://celestrak.com/NORAD/elements/weather.txt',
            'http://celestrak.com/NORAD/elements/resource.txt',
            'https://www.celestrak.com/NORAD/elements/cubesat.txt',
            'http://celestrak.com/NORAD/elements/stations.txt',
            'https://www.celestrak.com/NORAD/elements/sarsat.txt',
            'https://www.celestrak.com/NORAD/elements/noaa.txt',
            'https://www.celestrak.com/NORAD/elements/amateur.txt',
            'https://www.celestrak.com/NORAD/elements/engineering.txt')


LOGGER = logging.getLogger(__name__)
PKG_CONFIG_DIR = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'etc')


def read_platform_numbers(in_upper=False, num_as_int=False):
    """Read platform numbers from $PPP_CONFIG_DIR/platforms.txt."""
    out_dict = {}
    os.getenv('PPP_CONFIG_DIR', PKG_CONFIG_DIR)
    platform_file = None
    if 'PPP_CONFIG_DIR' in os.environ:
        platform_file = os.path.join(os.environ['PPP_CONFIG_DIR'], 'platforms.txt')
    if not platform_file or not os.path.isfile(platform_file):
        platform_file = os.path.join(PKG_CONFIG_DIR, 'platforms.txt')

    try:
        fid = open(platform_file, 'r')
    except IOError:
        LOGGER.error("Platform file %s not found.", platform_file)
        return out_dict
    for row in fid:
        # skip comment lines
        if not row.startswith('#'):
            parts = row.split()
            if len(parts) < 2:
                continue
            if in_upper:
                parts[0] = parts[0].upper()
            if num_as_int:
                parts[1] = int(parts[1])
            out_dict[parts[0]] = parts[1]
    fid.close()

    return out_dict


SATELLITES = read_platform_numbers(in_upper=True, num_as_int=False)
"""
The platform numbers are given in a file $PPP_CONFIG/platforms.txt
in the following format:

.. literalinclude:: ../../etc/platforms.txt
  :language: text
  :lines: 4-
"""


def read(platform, tle_file=None, line1=None, line2=None):
    """Read TLE for *platform*.

    The data are read from *tle_file*, from *line1* and *line2*, from
    the newest file provided in the TLES pattern, or from internet if
    none is provided.

    """
    return Tle(platform, tle_file=tle_file, line1=line1, line2=line2)


def get_all_platforms(tle_file=None):
    """Read all *platform* name from tle file
    tle_file - list path files
    """
    platforms = list()
    fid = io.open(tle_file, 'rb')
    for line in fid:
        line = line.decode('utf-8')
        line = line.strip()
        if not line.startswith('1') and not line.startswith('2') and len(line) > 0:
            platforms.append(line)
            l_1 = next(fid).decode('utf-8')
            l_2 = next(fid).decode('utf-8')
    fid.close()

    return platforms


def get_platform(tle_file=None):
    """Read *platform* name from tle file
    tle_file - path to tle file
    """
    try:
        fid = io.open(tle_file, 'rb')
        for line in fid:
            line = line.decode('utf-8')
            line = line.strip()
            if not line.startswith('1') and not line.startswith('2') and len(line) > 0:
                yield line
                l_1 = next(fid).decode('utf-8')
                l_2 = next(fid).decode('utf-8')
    except IOError:
        LOGGER.error("TLE file %s not found.", tle_file)
    finally:
        if fid:
            fid.close()


def fetch(destination):
    """Fetch TLE from internet and save it to `destination`."""
    with io.open(destination, mode="w", encoding="utf-8") as dest:
        for url in TLE_URLS:
            response = urlopen(url)
            dest.write(response.read().decode("utf-8"))


class ChecksumError(Exception):
    """ChecksumError."""

    pass


class Tle(object):
    """Class holding TLE objects."""

    def __init__(self, platform, tle_file=None, line1=None, line2=None):
        """Init."""
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
        """Return first TLE line."""
        return self._line1

    @property
    def line2(self):
        """Return second TLE line."""
        return self._line2

    @property
    def platform(self):
        """Return satellite platform name."""
        return self._platform

    def _checksum(self):
        """Calculate checksum for the current TLE."""
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
        """Read TLE data."""
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
                    if l_0.strip().upper() == self._platform:
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
        """Parse values from TLE data."""
        def _read_tle_decimal(rep):
            """Convert *rep* to decimal value."""
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
        self.epoch = \
            np.datetime64(dt.datetime.strptime(self.epoch_year, "%y") +
                          dt.timedelta(days=self.epoch_day - 1), 'us')
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
        """Format the class data for printing."""
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


PLATFORM_NAMES_TABLE = "(satid text primary key, platform_name text)"
SATID_TABLE = ("'{}' (epoch date primary key, tle text, insertion_time date,"
               " source text)")
SATID_VALUES = "INSERT INTO '{}' VALUES (?, ?, ?, ?)"
PLATFORM_VALUES = "INSERT INTO platform_names VALUES (?, ?)"
ISO_TIME_FORMAT = "%Y-%m-%dT%H:%M:%S.%f"


class Downloader(object):
    """Class for downloading TLE data."""

    def __init__(self, config):
        """Init."""
        self.config = config

    def fetch_plain_tle(self):
        """Fetch plain text-formated TLE data."""
        tles = {}
        if "fetch_plain_tle" in self.config["downloaders"]:
            sources = self.config["downloaders"]["fetch_plain_tle"]
            for source in sources:
                tles[source] = []
                failures = []
                for uri in sources[source]:
                    req = requests.get(uri)
                    if req.status_code == 200:
                        tles[source] += self.parse_tles(req.text)
                    else:
                        failures.append(uri)
                if len(failures) > 0:
                    logging.error(
                        "Could not fetch TLEs from %s, %d failure(s): [%s]",
                        source, len(failures), ', '.join(failures))
                logging.info("Downloaded %d TLEs from %s",
                             len(tles[source]), source)
        return tles

    def fetch_spacetrack(self):
        """Fetch TLE data from Space-Track."""
        tles = []
        login_url = "https://www.space-track.org/ajaxauth/login"
        download_url = ("https://www.space-track.org/basicspacedata/query/"
                        "class/tle_latest/ORDINAL/1/NORAD_CAT_ID/%s/format/"
                        "tle")
        download_url = download_url % ','.join(
            [str(key) for key in self.config['platforms']])

        user = self.config["downloaders"]["fetch_spacetrack"]["user"]
        password = self.config["downloaders"]["fetch_spacetrack"]["password"]
        credentials = {"identity": user, "password": password}

        with requests.Session() as session:
            # Login
            req = session.post(login_url, data=credentials)

            if req.status_code != 200:
                logging.error("Could not login to Space-Track")
                return tles

            # Get the data
            req = session.get(download_url)

            if req.status_code == 200:
                tles += self.parse_tles(req.text)
            else:
                logging.error("Could not retrieve TLEs from Space-Track")

        logging.info("Downloaded %d TLEs from %s", len(tles), "spacetrack")

        return tles

    def read_tle_files(self):
        """Read TLE data from files."""
        paths = self.config["downloaders"]["read_tle_files"]["paths"]

        # Collect filenames
        fnames = []
        for path in paths:
            if '*' in path:
                fnames += glob.glob(path)
            else:
                if not os.path.exists(path):
                    logging.error("File %s doesn't exist.", path)
                    continue
                fnames += [path]

        tles = []
        for fname in fnames:
            with open(fname, 'r') as fid:
                data = fid.read()
            tles += self.parse_tles(data)

        logging.info("Loaded %d TLEs from local files", len(tles))

        return tles

    def parse_tles(self, raw_data):
        """Parse all the TLEs in the given raw text data."""
        tles = []
        line1, line2 = None, None
        raw_data = raw_data.split('\n')
        for row in raw_data:
            if row.startswith('1 '):
                line1 = row
            elif row.startswith('2 '):
                line2 = row
            else:
                continue
            if line1 is not None and line2 is not None:
                try:
                    tle = Tle('', line1=line1, line2=line2)
                except ValueError:
                    logging.warning(
                        "Invalid data found - line1: %s, line2: %s",
                        line1, line2)
                else:
                    tles.append(tle)
                line1, line2 = None, None
        return tles


class SQLiteTLE(object):
    """Store TLE data in a sqlite3 database."""

    def __init__(self, db_location, platforms, writer_config):
        """Init."""
        self.db = sqlite3.connect(db_location)
        self.platforms = platforms
        self.writer_config = writer_config
        self.updated = False

        # Create platform_names table if it doesn't exist
        if not table_exists(self.db, "platform_names"):
            cmd = "CREATE TABLE platform_names " + PLATFORM_NAMES_TABLE
            with self.db:
                self.db.execute(cmd)
                logging.info("Created database table 'platform_names'")

    def update_db(self, tle, source):
        """Update the collected data.

        Only data with newer epoch than the existing one is used.

        """
        num = int(tle.satnumber)
        if num not in self.platforms:
            return
        tle.platform_name = self.platforms[num]
        if not table_exists(self.db, num):
            cmd = "CREATE TABLE " + SATID_TABLE.format(num)
            with self.db:
                self.db.execute(cmd)
                logging.info("Created database table '%d'", num)
            cmd = ""
            with self.db:
                self.db.execute(PLATFORM_VALUES, (num, self.platforms[num]))
                logging.info("Added platform name '%s' for ID '%d'",
                             self.platforms[num], num)
        cmd = SATID_VALUES.format(num)
        epoch = tle.epoch.item().isoformat()
        tle = '\n'.join([tle.line1, tle.line2])
        now = dt.datetime.utcnow().isoformat()
        try:
            with self.db:
                self.db.execute(cmd, (epoch, tle, now, source))
                logging.info("Added TLE for %d (%s), epoch: %s, source: %s",
                             num, self.platforms[num], epoch, source)
                self.updated = True
        except sqlite3.IntegrityError:
            pass

    def write_tle_txt(self):
        """Write TLE data to a text file."""
        if not self.updated and not self.writer_config.get('write_always',
                                                           False):
            return
        pattern = os.path.join(self.writer_config["output_dir"],
                               self.writer_config["filename_pattern"])
        now = dt.datetime.utcnow()
        fname = now.strftime(pattern)
        out_dir = os.path.dirname(fname)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            logging.info("Created directory %s", out_dir)
        data = []

        for satid, platform_name in self.platforms.items():
            if self.writer_config.get("write_name", False):
                data.append(platform_name)
            query = ("SELECT epoch, tle FROM '%s' ORDER BY "
                     "epoch DESC LIMIT 1" % satid)
            epoch, tle = self.db.execute(query).fetchone()
            date_epoch = dt.datetime.strptime(epoch, ISO_TIME_FORMAT)
            tle_age = (
                dt.datetime.utcnow() - date_epoch).total_seconds() / 3600.
            logging.info("Latest TLE for '%s' (%s) is %d hours old.",
                         satid, platform_name, int(tle_age))
            data.append(tle)

        with open(fname, 'w') as fid:
            fid.write('\n'.join(data))

        logging.info("Wrote %d TLEs to %s", len(data), fname)

    def close(self):
        """Close the database."""
        self.db.close()


def table_exists(db, name):
    """Check if the table 'name' exists in the database."""
    name = str(name)
    query = "SELECT 1 FROM sqlite_master WHERE type='table' and name=?"
    return db.execute(query, (name,)).fetchone() is not None


def main():
    """Run a test TLE reading."""
    tle_data = read('Noaa-19')
    print(tle_data)


if __name__ == '__main__':
    main()
