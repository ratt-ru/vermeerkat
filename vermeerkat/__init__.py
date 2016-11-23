#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2016 SKA South Africa
#
# This file is part of VerMeerKAT.
#
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

import logging
import os

import vermeerkat

from version import __version__

# Where is the module installed?
__install_path = os.path.split(os.path.abspath(vermeerkat.__file__))[0]

def create_logger():
    """ Create a console logger """
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    cfmt = logging.Formatter(('%(name)s - %(levelname)s - %(message)s'))

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(cfmt)

    log.addHandler(console)

    return log

# Create the log object
log = create_logger()

def install_path():
    """ Return the install path """
    return __install_path

def run(args):
    """ Execute the program """

    from config import configuration
    from observation import (query_recent_observations,
        download_observation)

    cfg = configuration(args)

    observations = query_recent_observations(cfg.solr_url)

    for o in observations:
        download_observation(o)