#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011 SMHI

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

setup(name='pyorbital',
      version="0.2.0",
      description='Orbital parameters and astronomical computations in Python',
      author='Martin Raspaud, Esben S. Nielsen',
      author_email='martin.raspaud@smhi.se, esn@dmi.dk',
      package_dir = {'pyorbital': 'pyorbital'},
      packages = ['pyorbital'],      
      install_requires=['numpy'],
      zip_safe=False
      )



