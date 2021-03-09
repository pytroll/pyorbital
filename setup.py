#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2011-2021 Pytroll Community
#
# Author(s):
#
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Panu Lahtinen <panu.lahtinen@fmi.fi>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
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

"""Setup for pyorbital."""

import os
from setuptools import setup, find_packages
import versioneer

try:
    with open('./README.md', 'r') as fd:
        long_description = fd.read()
except IOError:
    long_description = ''


setup(name='pyorbital',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Orbital parameters and astronomical computations in Python',
      author='The Pytroll Team',
      author_email='pytroll@googlegroups.com',
      classifiers=["Development Status :: 5 - Production/Stable",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License v3 " +
                   "or later (GPLv3+)",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python",
                   "Topic :: Scientific/Engineering",
                   "Topic :: Scientific/Engineering :: Astronomy"],
      url="https://github.com/pytroll/pyorbital",
      long_description=long_description,
      long_description_content_type='text/markdown',
      test_suite='pyorbital.tests.suite',
      packages=find_packages(),
      package_data={'pyorbital': [os.path.join('etc', 'platforms.txt')]},
      scripts=['bin/fetch_tles.py', ],
      install_requires=['numpy>=1.19.0', 'scipy', 'requests'],
      zip_safe=False,
      )
