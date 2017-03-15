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
import math

def launch(cfg, INPUT, MSDIR, OUTPUT, **kwargs):
    plot_name = kwargs["plot_name"]
    targets = kwargs["targets"]
    target_fields = kwargs["target_fields"]
    recipe = stimela.Recipe("Phase amplitude 2GC Pipeline", ms_dir=MSDIR)

    # Initial selfcal loop
    for target_field, target in zip(target_fields, targets):
        # Extract sources in mfs clean image to build initial sky model
        imname_prefix = cfg.obs.basename + "_SC0_3_" + plot_name[target.name]
        imname_mfs = imname_prefix + "-MFS-image.fits"

        # Construct LSM from sky cleaned after first SC subtraction
        old_model_prefix = cfg.obs.basename + "_LSM0_" + plot_name[target.name]
        old_model_name = old_model_prefix + ".lsm.html"
        model_prefix = cfg.obs.basename + "_LSM1_" + plot_name[target.name]
        model_name = model_prefix + ".lsm.html"
        recipe.add("cab/pybdsm", "extract_sources_%d" % target_field,
                   {
                       "image": "%s:output" % imname_mfs,
                       "outfile": model_prefix + ".fits",
                       "thresh_pix": cfg.source_find1.thresh_pix,
                       "thresh_isl": cfg.source_find1.thresh_isl,
                       "port2tigger": True,
                       "clobber": True,
                   },
                   input=INPUT, output=OUTPUT,
                   label="source_find1_%d:: Extract sources from previous round of cal" % target_field)

        # Stitch wsclean channel images of the previous run into a cube
        cubename = cfg.obs.basename + "_SC0_" + plot_name[target.name] + "-CLEAN_cube.fits"
        recipe.add("cab/fitstool", "fitstool",
                   {
                       "image": ['%s-%04d-image.fits:output' % (imname_prefix, a) for a in range(cfg.obs.im_numchans)],
                       "output": cubename,
                       "stack": True,
                       "fits-axis": 3,
                   },
                   input=INPUT, output=OUTPUT,
                   label="stitch_cube1_%d:: Stitch MFS image slices into a cube" % target_field)

        # Add SPIs
        recipe.add("cab/specfit", "add_SPIs_LSM0",
                   {
                       "image": "%s:output" % cubename,
                       "output-spi-image": "%s-spi.fits" % imname_prefix,
                       "output-spi-error-image": "%s-spi.error.fits" % imname_prefix,
                       "input-skymodel": "%s:output" % model_name,  # model to which SPIs must be added
                       "output-skymodel": model_name,
                       "tolerance-range": cfg.specfit1.tol,
                       "freq0": cfg.obs.freq_0,  # reference frequency for SPI calculation
                       "sigma-level": cfg.specfit1.sigma,
                   },
                   input=INPUT, output=OUTPUT,
                   label="SPI1_%d::Add SPIs to LSM" % target_field)

        # Selfcal and subtract brightest sources
        recipe.add("cab/calibrator", "Initial_Gjones_subtract_LSM1",
                   {
                       "skymodel": "%s:output" % model_name,
                       "label": cfg.selfcal1.label,
                       "msname": cfg.obs.msfile,
                       "threads": cfg.selfcal1.ncpu,
                       "column": cfg.selfcal1.column,
                       "output-data": cfg.selfcal1.output,
                       "Gjones": cfg.selfcal1.gjones,
                       "Gjones-solution-intervals": [cfg.selfcal1.gjones_time_interval,
                                                     cfg.selfcal1.gjones_freq_interval],
                       "DDjones-smoothing-intervals": cfg.selfcal1.ddjones_smoothing,
                       # TODO: MeerKAT beams need to go in this section
                       "Ejones": cfg.selfcal1.ejones,
                       "beam-files-pattern": cfg.selfcal1.beam_files_pattern,
                       "beam-l-axis": cfg.selfcal1.beam_l_axis,
                       "beam-m-axis": cfg.selfcal1.beam_m_axis,
                       "Gjones-ampl-clipping": cfg.selfcal1.gjones_ampl_clipping,
                       "Gjones-ampl-clipping-low": cfg.selfcal1.gjones_ampl_clipping_low,
                       "Gjones-ampl-clipping-high": cfg.selfcal1.gjones_ampl_clipping_high,
                       "Gjones-matrix-type": cfg.selfcal1.gjones_matrix_type,
                       "make-plots": True,
                       "field-id": target_field,
                       "Gjones-thresh-sigma": cfg.selfcal1.gjones_thresh_sigma,
                       "Gjones-chisq-clipping": cfg.selfcal1.gjones_chisq_clipping,
                       "tile-size" : cfg.selfcal1.tile_size
                   },
                   input=INPUT, output=OUTPUT,
                   label="SELFCAL1_%d:: Calibrate and subtract LSM0" % target_field)
        # Move the current flags to a new flag set in prep for next selfcal
        # round
        # recipe.add("cab/flagms", "backup_sc0_flags",
        #     {
        #         "msname"        : cfg.obs.msfile,
        #         "flagged-any"   : cfg.flagset_saveas_legacy.flagged_any,
        #         "flag"          : cfg.flagset_saveas_legacy.flag,
        #     },
        #     input=INPUT,   output=OUTPUT,
        #     label="flagset_selfcal0_flags:: Backup selfcal flags")

    # Initial selfcal loop
    phaseamp_2gc = []
    for target_field in target_fields:
        phaseamp_2gc += ["source_find1_%d" % target_field,
                         "stitch_cube1_%d" % target_field,
                         "SPI1_%d" % target_field,
                         "stitch_lsms1_%d" % target_field,
                         "SELFCAL1_%d" % target_field,
                         ]
    recipe.run(phaseamp_2gc)
