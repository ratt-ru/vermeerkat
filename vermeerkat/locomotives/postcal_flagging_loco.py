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
    bpcal_field = kwargs["bpcal_field"]
    gaincal_field = kwargs["gaincal_field"]
    target_fields = kwargs["target_fields"]
    recipe = stimela.Recipe("Post 1.1GC Flagging Engine", ms_dir=MSDIR)
    # Lets try squash the last of the RFI with the autoflagger
    # TODO: the strategy used here may still need some fine tuning
    # Think we should make it a bit more aggressive
    recipe.add("cab/autoflagger", "auto_flag_rfi_corrected_vis",
        {
            "msname"    : cfg.obs.msfile,
            "column"    : cfg.autoflag_corrected_vis.column,
            "strategy"  : cfg.autoflag_corrected_vis.strategy_file,
        },
        input=INPUT, output=OUTPUT,
        label="autoflag_corrected_vis:: Auto Flagging calibrated visibilities")

    # Some antennae may malfunction during the observation and not track
    # properly. This causes severe phase problems so its better to remove
    # them from the equation... for this we look at the phase of the calibrated
    # gain calibrator and flag out baselines (per channel) which are not up
    # to spec
    recipe.add("cab/politsiyakat", "flag_phases_amplitudes",
        {
            "task"                   : cfg.flag_phases_amplitudes.task,
            "msname"                 : cfg.obs.msfile,
            "data_column"            : cfg.flag_phases_amplitudes.data_column,
            "field"                  : ",".join(str(f) for f in target_fields +
                                                                [bpcal_field,
                                                                 gaincal_field]),
            "cal_field"              : ",".join(str(f) for f in [bpcal_field,
                                                                 gaincal_field]),
            "phase_range_clip"       : cfg.flag_phases_amplitudes.valid_phase_range,
            "invalid_count_frac_clip" : cfg.flag_phases_amplitudes.max_invalid_datapoints,
            "amp_frac_clip"          : cfg.flag_phases_amplitudes.amp_frac_clip,
            "output_dir"             : "",
            "nrows_chunk"            : cfg.flag_phases_amplitudes.nrows_chunk,
            "simulate"               : cfg.flag_phases_amplitudes.simulate,
            "nthreads"               : cfg.flag_phases_amplitudes.nthreads,
        },
        input=INPUT, output=OUTPUT,
        label="flag_phases_amplitudes:: Flag baselines based on calibrator phases")

    post_1gc_flagging = [ "autoflag_corrected_vis",
                          #"flag_phases_amplitudes"
                        ]
    recipe.run(post_1gc_flagging)
