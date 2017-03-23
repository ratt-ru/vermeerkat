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

def launch(cfg, INPUT, MSDIR, OUTPUT, **kwargs):
    recipe = stimela.Recipe("Conversion Engine", ms_dir=MSDIR)
    # Convert
    recipe.add("cab/h5toms", "h5toms",
        {
            'hdf5files'  : [cfg.obs.h5file],
            'output-ms'  : cfg.obs.msfile,
            'model-data' : True,
            'full_pol'   : cfg.h5toms.full_pol,
            'pols-to-use': cfg.h5toms.pols_to_use
        },
        input=INPUT, output=OUTPUT,
        label="convert::h5toms")

    # @mauch points out some MeerKAT h5 files contain flipped uvw
    # coordinates. Better to recalculate this from ANTENNAS
    # and be dead sure things will always work!
    recipe.add("cab/casa_fixvis", "fixvis",
        {
            "vis"       :   cfg.obs.msfile,
            "outputvis" :   cfg.obs.msfile, # into same ms please
            "reuse"     :   cfg.casa_fixvis.reuse,
        },
        input=INPUT, output=OUTPUT,
        label="recompute_uvw:: Recompute MeerKAT uvw coordinates")
    conversion_pipe = ["convert",
                       "recompute_uvw"]

    recipe.run(conversion_pipe)
