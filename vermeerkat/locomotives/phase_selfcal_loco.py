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
    recipe = stimela.Recipe("Initial phase-only 2GC Pipeline", ms_dir=MSDIR)

    # Add bitflag column. To keep track of flagsets
    recipe.add("cab/msutils", "msutils",
               {
                   'command': 'prep',
                   'msname': cfg.obs.msfile,
               },
               input=INPUT, output=OUTPUT,
               label="prepms::Adds flagsets")

    # Copy CORRECTED_DATA to DATA, so we can start selfcal
    recipe.add("cab/msutils", "shift_columns",
               {
                   "command": cfg.move_corrdata_to_data.command,
                   "msname": cfg.obs.msfile,
                   "fromcol": cfg.move_corrdata_to_data.fromcol,
                   "tocol": cfg.move_corrdata_to_data.tocol,
               },
               input=INPUT, output=OUTPUT,
               label="move_corrdata_to_data::msutils")

    # Initial selfcal loop
    for target_field, target in zip(target_fields, targets):
        # Extract sources in mfs clean image to build initial sky model
        imname_prefix = cfg.obs.basename + "_1GC_" + plot_name[target.name]
        imname_mfs = imname_prefix + "-MFS-image.fits"
        model_prefix = cfg.obs.basename + "_LSM0_" + plot_name[target.name]
        model_name = model_prefix + ".lsm.html"
        recipe.add("cab/pybdsm", "extract_sources_%d" % target_field,
                   {
                       "image": "%s:output" % imname_mfs,
                       "outfile": model_prefix + ".fits",
                       "thresh_pix": cfg.source_find0.thresh_pix,
                       "thresh_isl": cfg.source_find0.thresh_isl,
                       "port2tigger": True,
                   },
                   input=INPUT, output=OUTPUT,
                   label="source_find0_%d:: Extract sources from previous round of cal" % target_field)

        # Stitch wsclean channel images into a cube
        cubename = cfg.obs.basename + "_1GC_" + plot_name[target.name] + "-CLEAN_cube.fits"
        recipe.add("cab/fitstool", "fitstool",
                   {
                       "image": ['%s-%04d-image.fits:output' % (imname_prefix, a) for a in range(cfg.obs.im_numchans)],
                       "output": cubename,
                       "stack": True,
                       "fits-axis": 3,
                   },
                   input=INPUT, output=OUTPUT,
                   label="stitch_cube0_%d:: Stitch MFS image slices into a cube" % target_field)

        # Add SPIs
        recipe.add("cab/specfit", "add_SPIs_LSM0",
                   {
                       "image": "%s:output" % cubename,
                       "output-spi-image": "%s-spi.fits" % imname_prefix,
                       "output-spi-error-image": "%s-spi.error.fits" % imname_prefix,
                       "input-skymodel": "%s:output" % model_name,  # model to which SPIs must be added
                       "output-skymodel": model_name,
                       "tolerance-range": cfg.specfit0.tol,
                       "freq0": cfg.obs.freq_0,  # reference frequency for SPI calculation
                       "sigma-level": cfg.specfit0.sigma,
                   },
                   input=INPUT, output=OUTPUT,
                   label="SPI0_%d::Add SPIs to LSM" % target_field)

        # Selfcal and subtract brightest sources
        recipe.add("cab/calibrator", "Initial_Gjones_subtract_LSM0",
                   {
                       "skymodel": "%s:output" % model_name,
                       "label": cfg.selfcal0.label,
                       "msname": cfg.obs.msfile,
                       "threads": cfg.selfcal0.ncpu,
                       "column": cfg.selfcal0.column,
                       "output-data": cfg.selfcal0.output,
                       "Gjones": cfg.selfcal0.gjones,
                       "Gjones-solution-intervals": [int(math.ceil(float(cfg.obs.gain_sol_int[:-1]))),
                                                     int(cfg.obs.nchans /
                                                         float(100))],
                       "DDjones-smoothing-intervals": cfg.selfcal0.ddjones_smoothing,
                       # TODO: MeerKAT beams need to go in this section
                       "Ejones": cfg.selfcal0.ejones,
                       "beam-files-pattern": cfg.selfcal0.beam_files_pattern,
                       "beam-l-axis": cfg.selfcal0.beam_l_axis,
                       "beam-m-axis": cfg.selfcal0.beam_m_axis,
                       "Gjones-ampl-clipping": cfg.selfcal0.gjones_ampl_clipping,
                       "Gjones-ampl-clipping-low": cfg.selfcal0.gjones_ampl_clipping_low,
                       "Gjones-ampl-clipping-high": cfg.selfcal0.gjones_ampl_clipping_high,
                       "Gjones-matrix-type": cfg.selfcal0.gjones_matrix_type,
                       "make-plots": True,
                       "field-id": target_field
                   },
                   input=INPUT, output=OUTPUT,
                   label="SELFCAL0_%d:: Calibrate and subtract LSM0" % target_field)

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

        # create a mask for this round of selfcal
        maskname = imname_prefix + "_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask",
                   {
                       "image": '%s:output' % imname_mfs,
                       "output": maskname,
                       "sigma": cfg.cleanmask0.sigma,
                       "iters": cfg.cleanmask0.iters,
                       "boxes": cfg.cleanmask0.kernel,
                   },
                   input=INPUT, output=OUTPUT,
                   label="MSK_SC0_%d::Make clean mask" % target_field)

        # make another mfs image
        imname_prefix = cfg.obs.basename + "_SC0_" + plot_name[target.name]
        imname_mfs = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC0_%d" % target_field,
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_selfcal0.column,
                       "weight": "briggs %2.f" % (cfg.wsclean_selfcal0.robust),
                       "npix": cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_selfcal0.clean_iterations,
                       "mgain": cfg.wsclean_selfcal0.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_selfcal0.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "fitsmask": '%s:output' % maskname,
                       "auto-threshold": cfg.wsclean_selfcal0.autothreshold,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_SC0_%d::wsclean" % target_field)

        # full restore into the mfs
        imname_mfs_fullrest = imname_prefix + "-MFS-image.full_restore.fits"
        recipe.add("cab/tigger_restore", "full_restore_SC0_%d" % target_field,
                   {
                       "input-image": "%s:output" % imname_mfs,
                       "input-skymodel": "%s:output" % model_name,
                       "output-image": '%s:output' % imname_mfs_fullrest,
                       "freq": cfg.obs.freq_0,
                   },
                   input=INPUT, output=OUTPUT,
                   label="fullrestore_SC0_%d::tigger_restore" % target_field)

    # Initial selfcal loop
    phase_2gc = ["prepms",
                 "move_corrdata_to_data",
                 ]
    for target_field in target_fields:
        phase_2gc += ["source_find0_%d" % target_field,
                      "stitch_cube0_%d" % target_field,
                      "SPI0_%d" % target_field,
                      "SELFCAL0_%d" % target_field,
                      "MSK_SC0_%d" % target_field,
                      "image_SC0_%d" % target_field,
                      "fullrestore_SC0_%d" % target_field,
                      ]

    recipe.run(phase_2gc)