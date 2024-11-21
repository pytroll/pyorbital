# -*- coding: utf-8 -*-

# Copyright (c) 2017-2024 Pytroll Community

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

"""Package file."""

import numpy as np

from pyorbital.version import __version__  # noqa


def dt2np(utc_time):
    """Convert datetime to numpy datetime64 object."""
    try:
        return np.datetime64(utc_time)
    except ValueError:
        return utc_time.astype("datetime64[ns]")
