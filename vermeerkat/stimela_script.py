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

import ast
import os
import sys

import stimela
import vermeerkat

from vermeerkat.config import configuration
from vermeerkat.observation import (query_recent_observations,
    download_observation)

INPUT = "input"
OUTPUT = "output"
PREFIX = "vermeerkat-pipeline"
MSDIR  = "msdir"

# So that we can access GLOBALS pass through to the run command
stimela.register_globals()

# args is a global string variable representing a python list of strings
# Evaluate it to get the list of strings back
args = ast.literal_eval(args)

# Load in the configuration
cfg = configuration(args)

# If a specific HDF5 file is specified, use that to specify
# the observation query
if cfg.general.hdf5_file is not None:
    query = 'Filename:{}'.format(cfg.general.hdf5_file)
else:
    query = None

# Find recent observations
observations = query_recent_observations(cfg.general.solr_url, query)

if len(observations) == 0:
    vermeerkat.log.warn("No observations matching query string '{}'".format(query))

# Image each observation
for o in observations:
    f = download_observation(o, INPUT)
    h5file = os.path.split(f)[1]
    basename = os.path.splitext(h5file)[0]
    msfile = ''.join((basename, '.ms'))
    msname = ''.join((basename, '_1284.full_pol.ms'))

    recipe = stimela.Recipe("Imaging Pipeline", ms_dir=MSDIR)
    recipe.add("cab/h5toms", "h5toms",
        {
            'hdf5files' : [h5file],
            'output-ms' : msfile
        },
        input=INPUT, output=OUTPUT,
        label="convert::h5toms")

    recipe.add("cab/rfimasker", "mask_stuff",
        {
            "msname" : msname,
            "mask"   : 'rfi_mask.pickle',
        },
        input=INPUT, output=OUTPUT,
        label="mask::maskms")

    recipe.add("cab/autoflagger", "auto_flag_rfi",
        {
            "msname"    : msname,
            "column"    : "DATA",
            #"strategy"  : <RFI strategy file>,
        },
        input=INPUT, output=OUTPUT,
        label="autoflag1:: Auto Flagging ms")

    recipe.add("cab/wsclean", "wsclean",
        {
            "msname"            : msname,
            "column"            : 'DATA',
            "weight"            : 'briggs',
            "robust"            : 0,
            "npix"              : 4096,
            "cellsize"          : 1,
            "clean_iterations"  : 1000,
            "mgain"             : 0.9,
        },
        input=INPUT, output=OUTPUT,
        label="image::wsclean")

    recipe.run()
