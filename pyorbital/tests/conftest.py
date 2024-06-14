#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Pyorbital developers

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

"""Fixtures and helper functions for unittests."""

import pytest


TEST_YAML_CONFIG_CONTENT = """
station:
  longitude: 16.15
  latitude: 58.58
  altitude: 0.0

# Here a list of directories can be given, to determine where to look up the
# TLE files covering the time period of interest.
tle-dirs:
  - /path/to/tle/files

# The format of tle filenames.
tle-file-format: 'TLE_{platform:s}.txt'
"""


@pytest.fixture
def fake_yamlconfig_file(tmp_path):
    """Write fake yaml config file."""
    file_path = tmp_path / 'snos_mystation.yaml'
    with open(file_path, 'w') as fpt:
        fpt.write(TEST_YAML_CONFIG_CONTENT)

    yield file_path


TEST_TLE_FILE1_CALIPSO = """1 29108U 06016B   13365.56860824  .00000817  00000-0  19114-3 0  5182
2 29108  98.2183 306.3131 0001262  81.5858 278.5464 14.57147103408355
1 29108U 06016B   14001.73596105  .00000831  00000-0  19433-3 0  5198
2 29108  98.2189 307.4643 0001241  81.3032 278.8307 14.57149241408523
1 29108U 06016B   14002.83464415  .00000945  00000-0  21969-3 0  5203
2 29108  98.2193 308.5492 0001225  81.5021 278.6342 14.57151825408689
1 29108U 06016B   14004.00199336 -.00003567  00000-0 -78180-3 0  5214
2 29108  98.2205 309.7016 0001599  74.2235 285.7542 14.57124395408850
"""

TEST_TLE_FILE1_SNPP = """1 37849U 11061A   13365.52212068  .00000110  00000-0  73158-4 0  6643
2 37849  98.7742 301.9172 0001322 129.8543 353.2962 14.19526864112801
1 37849U 11061A   14001.84985892  .00000068  00000-0  53438-4 0  6650
2 37849  98.7742 303.2326 0001351 129.0203 295.4516 14.19527275112993
1 37849U 11061A   14002.47164083  .00000106  00000-0  50370-4 0  6667
2 37849  98.7717 303.8433 0001170 118.9961 241.1812 14.19527267113086
1 37849U 11061A   14003.12858987  .00000112  00000-0  74213-4 0  6662
2 37849  98.7741 304.4994 0001384 128.2427 347.2535 14.19527795113176
1 37849U 11061A   14004.18951829  .00000140  00000-0  87286-4 0  6674
2 37849  98.7741 305.5505 0001426 126.5558   7.5558 14.19528199113326
"""


@pytest.fixture
def fake_tle_file1_calipso(tmp_path):
    """Write fake TLE file for Calipso."""
    file_path = tmp_path / 'TLE_calipso.txt'
    with open(file_path, 'w') as fpt:
        fpt.write(TEST_TLE_FILE1_CALIPSO)

    yield file_path


@pytest.fixture
def fake_tle_file1_snpp(tmp_path):
    """Write fake TLE file for Suomi-NPP."""
    file_path = tmp_path / 'TLE_snpp.txt'
    with open(file_path, 'w') as fpt:
        fpt.write(TEST_TLE_FILE1_SNPP)

    yield file_path
