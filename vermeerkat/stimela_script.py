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

import stimela

INPUT = "input"
OUTPUT = "output"
PREFIX = "stimela-example"
MSDIR  = "msdir"

# So that we can access GLOBALS pass through to the run command
stimela.register_globals()

recipe = stimela.Recipe("Imaging Pipeline", ms_dir=MSDIR)
recipe.add("cab/h5toms",
    "h5toms",
    {'hdf5files' : [h5file] },
    input=INPUT,
    output=OUTPUT)

recipe.run()