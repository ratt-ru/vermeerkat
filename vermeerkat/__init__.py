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

__stimela_script_path = os.path.join(__install_path, 'stimela_script.py')

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

    observations = query_recent_observations(cfg.general.solr_url)

    def _create_stimela_dirs():
        # Create stimela input and output directories
        for path in ('input', 'output', 'msdir'):
            full_path = os.path.join(os.getcwd(), path)

            if not os.path.exists(full_path):
                os.mkdir(full_path)
            else:
                if not os.path.isdir(full_path):
                    raise ValueError("Error creating directory '{}'. "
                        "A file already exists at this location.")

            yield (path, full_path)

    import stimela

    stimela_dirs = { p: f for p, f in _create_stimela_dirs() }

    for o in observations:
        f = download_observation(o, stimela_dirs['input'])
        h5file = os.path.split(f)[1]
        basename = os.path.splitext(h5file)[0]
        msfile = ''.join((basename, '.ms'))

        stimela.run(['-g', 'h5file={}:str'.format(h5file),
            '-g', 'output_ms={}:str'.format(msfile),
            __stimela_script_path])
