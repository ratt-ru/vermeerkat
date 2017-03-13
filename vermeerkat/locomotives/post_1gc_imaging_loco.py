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
    targets = kwargs["targets"]
    field_index = kwargs["field_index"]
    target_fields = kwargs["target_fields"]
    recipe = stimela.Recipe("1GC Imaging Engine", ms_dir=MSDIR)
    # imaging
    for target_field, target in zip(target_fields, targets):
        imname = cfg.obs.basename + "_1GC_" + target.name
        recipe.add("cab/wsclean", "wsclean_%d" % field_index[target.name],
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image.column,
                       "weight": "briggs %.2f" % (cfg.wsclean_image.robust),
                       "npix": cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image.clean_iterations,
                       "mgain": cfg.wsclean_image.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image.joinchannels,
                       "field": target_field,
                       "name": imname,
                       "auto-threshold": cfg.wsclean_image.autothreshold,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_%d::wsclean" % target_field)

    recipe.run(["image_%d" % t for t in target_fields])