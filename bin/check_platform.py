# Copyright (c) 2023 Pyorbital Developers
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
"""Check if a satellite is supported on default.

If not the name and its NORAD number needs to be added to a local copy of the
platforms.txt file, which then needs to be placed in the directory pointed to
by the environment variable PYORBITAL_CONFIG_PATH.

"""

import argparse
import logging
from pyorbital.tlefile import check_is_platform_supported
from pyorbital.logger import logging_on

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Check if a satellite is supported.')
    parser.add_argument("-s", "--satellite",
                        help=("Name of the Satellite - following WMO Oscar naming."),
                        default=None,
                        required=True,
                        type=str)

    args = parser.parse_args()
    satellite_name = args.satellite

    logging_on(logging.INFO)
    check_is_platform_supported(satellite_name)
