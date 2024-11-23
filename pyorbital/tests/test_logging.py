#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2024 Pytroll Community

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

"""Test the logging module."""

import logging

from pyorbital.logger import get_logger, logging_off, logging_on


def test_logging_on_and_off(caplog):
    """Test that switching logging on and off works."""
    logger = get_logger("pyorbital.spam")
    logging_on()
    with caplog.at_level(logging.WARNING):
        logger.debug("I'd like to leave the army please, sir.")
        logger.warning("Stop that! It's SPAM.")
    assert "Stop that! It's SPAM" in caplog.text
    assert "I'd like to leave the army please, sir." not in caplog.text
    logging_off()
    with caplog.at_level(logging.DEBUG):
        logger.warning("You've got a nice army base here, Colonel.")
    assert "You've got a nice army base here, Colonel." not in caplog.text
