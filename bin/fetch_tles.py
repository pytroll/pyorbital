#!/usr/bin/env python

"""Script to download and store satellite TLE data."""

import sys
import os
import glob
import sqlite3
import datetime as dt
import logging
import logging.config

import yaml
import requests
from pyorbital.tlefile import Tle

PLATFORM_NAMES_TABLE = "(satid text primary key, platform_name text)"
SATID_TABLE = "'{}' (epoch date primary key, tle text, insertion_time date, source text)"
SATID_VALUES = "INSERT INTO '{}' VALUES (?, ?, ?, ?)"
PLATFORM_VALUES = "INSERT INTO platform_names VALUES (?, ?)"


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
                    logging.error("Could not fetch TLEs from %s, %d failure(s): [%s]",
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
                    continue
                fnames += path

        tles = []
        for fname in fnames:
            with open(fname, 'r') as fid:
                data = fid.read()
            tles += self.parse_tles(data)

        logging.info("Loaded %d TLEs from local files", len(tles))

        return tles

    def parse_tles(self, raw_data):
        """Parse all the TLEs in the given raw text data.

        Return only the platforms that are configured.

        """
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
                tle = Tle('', line1=line1, line2=line2)
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
        now = dt.datetime.utcnow()
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
        if not self.updated and not self.writer_config.get('write_always', False):
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
            query = "SELECT tle FROM '%s' ORDER BY epoch DESC LIMIT 1" % satid
            tle = self.db.execute(query).fetchone()[0]
            data.append(tle)

        with open(fname, 'w') as fid:
            fid.write('\n'.join(data))

        logging.info("Wrote %d TLEs to %s", len(data), fname)

    def close(self):
        """Close the database."""
        self.db.close()


def read_config(config_fname):
    """Read and parse config file."""
    with open(config_fname, 'r') as fid:
        config = yaml.load(fid, Loader=yaml.SafeLoader)
    return config


def table_exists(db, name):
    """Check if the table 'name' exists in the database."""
    name = str(name)
    query = "SELECT 1 FROM sqlite_master WHERE type='table' and name=?"
    return db.execute(query, (name,)).fetchone() is not None


def value_exists(db, table, name):
    """Check if the table 'name' exists in the database."""
    name = str(name)
    query = "SELECT 1 FROM sqlite_master WHERE type='table' and name=?"
    return db.execute(query, (name,)).fetchone() is not None


def main():
    """Run TLE downloader."""
    config = read_config(sys.argv[1])
    if 'logging' in config:
        logging.config.dictConfig(config['logging'])
    else:
        logging.basicConfig(level=logging.INFO)

    downloader = Downloader(config)
    db = SQLiteTLE(config['database']['path'], config['platforms'],
                   config['text_writer'])

    logging.info("Start downloading TLEs")
    for dl_ in config['downloaders']:
        fetcher = getattr(downloader, dl_)
        tles = fetcher()
        if isinstance(tles, dict):
            for source in tles:
                for tle in tles[source]:
                    db.update_db(tle, source)
        else:
            source = 'file'
            if "spacetrack" in dl_:
                source = 'spacetrack'
            for tle in tles:
                db.update_db(tle, source)

    db.write_tle_txt()
    db.close()
    logging.info("TLE downloading finished")


if __name__ == "__main__":
    main()
