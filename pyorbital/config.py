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

"""Helper functions to read and extract configuration parameters."""

import yaml
from yaml import SafeLoader
from collections.abc import Mapping


def _recursive_dict_update(d, u):
    """Recursive dictionary update.

    Copied from:

        http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth

    """
    for k, v in u.items():
        if isinstance(v, Mapping):
            r = _recursive_dict_update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d


def get_config(configfile):
    """Get the configuration from file.

    :configfile: The file path of the yaml configuration file.

    :return: A configuration dictionary.

    """
    config = {}
    with open(configfile, 'r') as fp_:
        config = _recursive_dict_update(config, yaml.load(fp_, Loader=SafeLoader))

    return config
