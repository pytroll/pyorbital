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

"""Test the SNO finder functionality."""


import pytest
import numpy as np
import datetime as dt
from datetime import timezone
from pyorbital.sno_utils import SNOfinder
from pyorbital.sno_utils import check_overlapping_times
# from pyorbital.config import get_config
from pyorbital.tests.test_helper import get_dataframe_from_ascii
from pyorbital.tests.test_helper import df2geojson

# from pyorbital.sno_utils import geojson_compare_derived_snos_against_reference
from pyorbital.sno_utils import get_closest_sno_to_reference

# The data below are until further regarded as the "truth".
# The SNOs here have been found by running a Fotran program provided by NOAA colleagues to SMHI around ~2008
# Adam Dybbroe, 2024-06-14
TEST1_SNOS_ASCII = """
  2014  1  3  4 41  40.9     2014  1  3  4 43  16.7        78.34    -0.56
  2014  1  3  5 32  15.5     2014  1  3  5 32  34.3       -78.74   169.19
  2014  1  3  6 22  51.0     2014  1  3  6 21  51.8        79.08   -21.21
"""

TEST2_SNOS_ASCII = """
  2014  1  3  4 41  40.9     2014  1  3  4 43  16.7        78.34    -0.56
  2014  1  3  5 32  15.5     2014  1  3  5 32  34.3       -78.74   169.19
  2014  1  3  6 22  51.0     2014  1  3  6 21  51.8        79.08   -21.21
  2014  1  5 20 58  34.8     2014  1  5 20 59  41.5        78.51   116.17
  2014  1  5 21 49  10.0     2014  1  5 21 48  59.9       -78.88   -74.14
  2014  1  5 22 39  45.8     2014  1  5 22 38  17.8        79.19    95.39
  2014  1  8 13 15  31.6     2014  1  8 13 16  25.4        78.55  -127.76
  2014  1  8 14  6   6.9     2014  1  8 14  5  43.8       -78.92    41.90
  2014  1  8 14 56  43.0     2014  1  8 14 55   2.2        79.23  -148.59
  2014  1 11  4 41  53.6     2014  1 11  4 43  42.6       -78.21   178.74
  2014  1 11  5 32  27.1     2014  1 11  5 32  58.4        78.64   -11.43
  2014  1 11  6 23   2.8     2014  1 11  6 22  17.1       -78.99   158.18
"""


def test_get_snos_calipso_snpp(fake_yamlconfig_file, fake_tle_file1_calipso, fake_tle_file1_snpp):
    """Test finding calipso and SNPP SNOs within time window."""
    platform_one = 'snpp'
    platform_two = 'calipso'
    time_window = (dt.datetime(2014, 1, 3, tzinfo=timezone.utc),
                   dt.datetime(2014, 1, 4, tzinfo=timezone.utc))
    sno_min_thr = 2  # SNOs allowed to have 2 minute deviation between the two platforms
    mysnofinder = SNOfinder(platform_one, platform_two, time_window, sno_min_thr)
    mysnofinder.set_configuration(fake_yamlconfig_file)
    mysnofinder._conf['tle-dirs'].insert(0, str(fake_tle_file1_snpp.parent))
    mysnofinder._conf['tle-dirs'].insert(0, str(fake_tle_file1_calipso.parent))

    results = mysnofinder.get_snos_within_time_window()

    np.testing.assert_array_equal(results['within_local_reception_area'].values,
                                  np.array([True, False, True]))

    mysnofinder.dataframe2geojson(results)

    test_ref = get_dataframe_from_ascii(TEST1_SNOS_ASCII)
    ref_features = df2geojson(test_ref)

    for idx, expected in enumerate([6.5, 12, 18]):
        eval_res = get_closest_sno_to_reference(mysnofinder.geojson_results, ref_features['features'][idx])
        assert eval_res[0] < expected
        assert check_overlapping_times(eval_res[1], eval_res[2])


# Times: 2014-01-03 06:22:53.516417 2014-01-03 06:21:54.322415
# ('2014-01-03 06:22:51', '2014-01-03 06:21:51.800000')
TWIN1 = (dt.datetime(2014, 1, 3, 6, 22, 51).replace(tzinfo=timezone.utc),
         dt.datetime(2014, 1, 3, 6, 21, 51, 800000).replace(tzinfo=timezone.utc))
TWIN2 = (dt.datetime(2014, 1, 3, 6, 22, 53, 516417).replace(tzinfo=timezone.utc),
         dt.datetime(2014, 1, 3, 6, 21, 51, 800000).replace(tzinfo=timezone.utc))


@pytest.mark.parametrize("twin1, twin2, expected",
                         [(TWIN1, TWIN2, True)
                          ]
                         )
def test_check_overlapping_times(twin1, twin2, expected):
    """Test the check to see if two time windows overlap."""
    result = check_overlapping_times(twin1, twin2)
    assert result
