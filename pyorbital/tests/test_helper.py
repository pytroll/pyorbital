#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c22526.ad.smhi.se>

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

"""Helper functions for unittesting."""

import pandas as pd
import io
from datetime import datetime, timedelta
# from geojson import dump


def get_dataframe_from_ascii(sno_ascii_data):
    """Get data as pandas dataframe from ascii test file."""
    # Define the column names
    column_names = [
        'year1', 'month1', 'day1', 'hour1', 'minute1', 'second1',
        'year2', 'month2', 'day2', 'hour2', 'minute2', 'second2',
        'latitude', 'longitude'
    ]
    # Define the widths of each column in the fixed-width format
    colspecs = [
        (0, 6), (6, 9), (9, 12), (12, 15), (15, 18), (18, 24),  # First datetime components
        (24, 33), (33, 36), (36, 39), (39, 42), (42, 45), (45, 51),  # Second datetime components
        (51, 64), (64, 73)  # Latitude and longitude
    ]

    str_data = io.StringIO(sno_ascii_data)
    # Read the data using pandas
    df = pd.read_fwf(str_data, colspecs=colspecs, header=None, names=column_names)
    return df


def df2geojson(df):
    """Convert dataframe to Geojson."""
    # geojson skeleton
    geojson = {"type": "FeatureCollection", "features": []}

    # go through dataframe, append entries to geojson format
    for _, row in df.iterrows():
        dtobj1 = datetime(int(row['year1']), int(row['month1']), int(row['day1']),
                          int(row['hour1']), int(row['minute1'])) + timedelta(seconds=float(row['second1']))
        dtobj2 = datetime(int(row['year2']), int(row['month2']), int(row['day2']),
                          int(row['hour2']), int(row['minute2'])) + timedelta(seconds=float(row['second2']))

        feature = {"type": "Feature", "geometry": {"type": "Point",
                                                   "coordinates": [row['longitude'],
                                                                   row['latitude']]},
                   "properties": {"datetime1": dtobj1.isoformat(),
                                  "datetime2": dtobj2.isoformat()}}
        geojson['features'].append(feature)

    # with open('test1_sno.geojson', 'w') as fp:
    #    json.dump(geojson, fp)

    return geojson
