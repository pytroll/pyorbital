#!/usr/bin/env python

"""Script to download and store satellite TLE data."""

import sys
import logging
import logging.config

import yaml
from pyorbital.tlefile import Downloader, SQLiteTLE


def read_config(config_fname):
    """Read and parse config file."""
    with open(config_fname, 'r') as fid:
        config = yaml.load(fid, Loader=yaml.SafeLoader)
    return config


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
