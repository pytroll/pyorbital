#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011-2014

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

from setuptools import setup
import imp

version = imp.load_source('pyorbital.version', 'pyorbital/version.py')

setup(name='pyorbital',
      version=version.__version__,
      description='Orbital parameters and astronomical computations in Python',
      author='Martin Raspaud, Esben S. Nielsen',
      author_email='martin.raspaud@smhi.se, esn@dmi.dk',
      classifiers=["Development Status :: 5 - Production/Stable",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License v3 " +
                   "or later (GPLv3+)",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python",
                   "Topic :: Scientific/Engineering",
                   "Topic :: Scientific/Engineering :: Astronomy"],
      url="https://github.com/mraspaud/pyorbital",
      test_suite='pyorbital.tests.suite',
      package_dir = {'pyorbital': 'pyorbital'},
      packages = ['pyorbital'],      
      install_requires=['numpy>=1.6.0'],
      zip_safe=False,
      )

