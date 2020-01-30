#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Martin Raspaud

# Author(s):

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
"""The tests package."""

from pyorbital.tests import (test_aiaa, test_tlefile, test_orbital,
                             test_astronomy, test_geoloc)
import unittest


def suite():
    """The global test suite."""
    mysuite = unittest.TestSuite()
    # Test the documentation strings
    # mysuite.addTests(doctest.DocTestSuite(image))
    # Use the unittests also
    mysuite.addTests(test_aiaa.suite())
    mysuite.addTests(test_tlefile.suite())
    mysuite.addTests(test_orbital.suite())
    mysuite.addTests(test_astronomy.suite())
    mysuite.addTests(test_geoloc.suite())
    return mysuite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
