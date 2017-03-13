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


import os
import copy
import math

import stimela
import vermeerkat
import vermeerkat.caltables as vmct
import vermeerkat.config as vmc
import vermeerkat.observation as vmo
import vermeerkat.utils as vmu
from vermeerkat.locomotives import converter_loco, \
                                   flagging_loco, \
                                   init_casa_1gc_loco, \
                                   postcal_flagging_loco, \
                                   secondpass_casa_1gc_loco, \
                                   post_1gc_diagnostics_loco

# So that we can access GLOBALS pass through to the run command
stimela.register_globals()

# Load in the configuration
cfg = vmc.configuration(vmc.retrieve_args())

# Register directories
INPUT = cfg.general.input
OUTPUT = cfg.general.output
PREFIX = cfg.general.prefix
MSDIR  = cfg.general.msdir

# Get list of custom calibrators
calibrator_db = vmct.calibrator_database()

# Get a list of observations
obs_metadatas = vmo.observation_metadatas(INPUT, cfg)

if len(obs_metadatas) == 0:
    vermeerkat.log.warn("No observations found for given parameters")

# Image each observation
for obs_metadata in obs_metadatas:
    # Merge observation metadata into our config
    vmu.merge_observation_metadata(cfg, obs_metadata)

    # Dump the observation metadata
    vmo.dump_observation_metadata(INPUT, ''.join([cfg.obs.basename, '.json']),
        obs_metadata)

    # Download if no valid observation exists
    if not vmo.valid_observation_exists(INPUT, cfg.obs.h5file, obs_metadata):
        vmo.download_observation(INPUT, obs_metadata)

    # Load in scans
    scans = vmu.load_scans(os.path.join(INPUT, cfg.obs.h5file))

    # Map scan target name to a list of scans associated with it
    field_scan_map = vmu.create_field_scan_map(scans)

    #Get the start and end of the observation
    start_obs, end_obs = vmu.start_end_observation(os.path.join(INPUT,
                                                                cfg.obs.h5file))

    # Categories the fields observed in each scan
    field_index, bpcals, gaincals, targets = vmu.categorise_fields(scans)

    # Log useful information
    vermeerkat.log.info("The following fields were observed:")

    # Sort fields by index
    for k, v in sorted(field_index.items(), key=lambda (k, v): v):
        tags = set.union(*(set(s.tags) for s in field_scan_map[k]))
        scan_seconds = sum(s.scan_length for s in field_scan_map[k])
        vermeerkat.log.info("\t %d: %s %s %s" % (v, k.ljust(20),
                    vmu.fmt_seconds(scan_seconds).ljust(20),
                    [t for t in tags]))

    # Use nicer names for source plots
    plot_name = { s: s.replace(' ', '_') for s
                                    in field_index.keys() }

    if len(targets) < 1:
        raise RuntimeError("Observation %s does not have any "
                           "targets" % obs_metadata["ProductName"])
    if len(gaincals) < 1:
        raise RuntimeError("Observation %s does not have any "
                           "gain calibrators" % obs_metadata["ProductName"])
    if len(bpcals) < 1:
        raise RuntimeError("Observation %s does not have any "
                           "bandpass calibrators" % obs_metadata["ProductName"])

    # Select a gain calibrator
    gaincal_index, gain_cal = vmu.select_gain_calibrator(cfg,
                                                         targets,
                                                         gaincals,
                                                         observation_start=start_obs,
                                                         scans=scans,
                                                         strategy=cfg.general.gcal_selection_heuristic)

    # Compute the observation time spent on gain calibrators
    gaincal_scan_times = vmu.total_scan_times(field_scan_map, gaincals)

    # Compute observation time on bandpass calibrators
    bandpass_scan_times = vmu.total_scan_times(field_scan_map, bpcals)
    # Select the bandpass calibrator
    bpcal_index, bandpass_cal = vmu.select_bandpass_calibrator(cfg,
                                        bpcals, bandpass_scan_times)

    # Choose the solution interval for the bandpass calibrator
    # by choosing the bandpass calibrator's minimum scan length
    bpcal_sol_int = min(s.scan_length for s in field_scan_map[bandpass_cal.name])

    # Compute observation time on target
    target_scan_times = vmu.total_scan_times(field_scan_map, targets)

    # Get field ids for bandpass, gaincal and target
    bpcal_field = field_index[bandpass_cal.name]
    gaincal_field = field_index[gain_cal.name]
    target_fields = [field_index[t.name] for t in targets]

    vermeerkat.log.info("The following targets were observed:")

    for target, target_scan_time in zip(targets, target_scan_times):
        vermeerkat.log.info("\t %s (%d) - total observation time %s." %
                                (target.name,
                                 field_index[target.name],
                                 vmu.fmt_seconds(scan_seconds)))

    vermeerkat.log.info("Using '%s' (%d) as a gain calibrator "
                        "(total observation time: %s)" %
                            (gain_cal.name,
                             field_index[gain_cal.name],
                             vmu.fmt_seconds(gaincal_scan_times[gaincal_index])))
    vermeerkat.log.info("Using '%s' (%d) as a bandpass calibrator (total "
                "observation time: %s, minimum scan time: %s)" %
                            (bandpass_cal.name,
                             field_index[bandpass_cal.name],
                             vmu.fmt_seconds(bandpass_scan_times[bpcal_index]),
                             vmu.fmt_seconds(bpcal_sol_int)))
    vermeerkat.log.info("Imaging the following targets: '%s' (%s)" %
                            (",".join([t.name for t in targets]),
                             ",".join([str(field_index[t.name]) for t in targets])))
    vermeerkat.log.info("Writing ms file to '%s'" % cfg.obs.msfile)
    vermeerkat.log.info("Using firstpass aoflagger strategy: %s" %
                            cfg.autoflag.strategy_file)
    vermeerkat.log.info("Using secondpass aoflagger strategy: %s" %
                            cfg.autoflag_corrected_vis.strategy_file)
    vermeerkat.log.info("Using RFI mask: %s" %
                            cfg.rfimask.rfi_mask_file)
    vermeerkat.log.info("Correlator integration interval "
                        "recorded as: %.2f secs" %
                            cfg.obs.correlator_integration_time)
    vermeerkat.log.info("MFS maps will contain %.3f MHz per slice" %
                            (cfg.obs.bw_per_image_slice / 1e6))
    vermeerkat.log.info("Maps will cover %.3f square degrees at "
                        "angular resolution %.3f asec" %
                            (cfg.obs.fov / 3600.0,
                            cfg.obs.angular_resolution))
    vermeerkat.log.warn("Assuming maximum baseline "
                        "is %.2f meters" %
                            cfg.obs.telescope_max_baseline)
    vermeerkat.log.info("Will use '%s' as reference antenna" %
                            cfg.obs.refant)
    vermeerkat.log.info("Observed band covers %.2f MHz "
                        "+/- %.2f MHz averaged into %d channels" %
                        (cfg.obs.freq_0 / 1e6,
                        (cfg.obs.nchans // 2) * cfg.obs.chan_bandwidth / 1e6,
                        cfg.obs.nchans))
    vermeerkat.log.info("Observation start time: %s UTC. Observation end time "
                        "%s UTC. Length %s" %
        (start_obs.strftime("%Y-%m-%d %H:%M:%S"),
         end_obs.strftime("%Y-%m-%d %H:%M:%S"),
         str(end_obs - start_obs)))

    #########################################################################
    #
    # Conversions and fixes
    #
    #########################################################################
    converter_loco.launch(cfg, INPUT, MSDIR, OUTPUT)

    #########################################################################
    #
    # Preliminary RFI, band and autocorr flagging
    #
    #########################################################################
    flagging_loco.launch(cfg, INPUT, MSDIR, OUTPUT)

    #########################################################################
    #
    # Initial 1GC
    #
    # We first calibrate the calibrators for flagging
    # purposes, we will want to discard these solutions when we
    # did mitigation flagging and generate, hopefully, improved 1GC
    # solutions
    #########################################################################
    init_casa_1gc_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                              bandpass_cal=bandpass_cal,
                              calibrator_db=calibrator_db,
                              bpcal_field=bpcal_field,
                              gaincal_field=gaincal_field,
                              bpcal_sol_int=bpcal_sol_int,
                              target_fields=target_fields)

    postcal_flagging_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                                 bpcal_field=bpcal_field,
                                 gaincal_field=gaincal_field,
                                 target_fields=target_fields)

    #########################################################################
    #
    # 1GC solutions gained from mitigation flagging
    #
    # Now that we mitigated some of the worst RFI and observation problems
    # we can generate better 1GC solutions
    #
    # This is an optional step if we cannot generate new solutions due
    # to significant flagging we can skip and use the corrected data
    # as is.
    #########################################################################
    secondpass_casa_1gc_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                                    bandpass_cal=bandpass_cal,
                                    calibrator_db=calibrator_db,
                                    bpcal_field=bpcal_field,
                                    gaincal_field=gaincal_field,
                                    bpcal_sol_int=bpcal_sol_int,
                                    target_fields=target_fields)

    #########################################################################
    #
    # Post 1GC diagnostic plots
    #
    # With our best solutions in hand we plot some results for the observer
    #########################################################################
    post_1gc_diagnostics_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                                     bandpass_cal=bandpass_cal,
                                     bpcal_field=bpcal_field,
                                     gaincal_field=gaincal_field)

    #########################################################################
    #
    # Post 1GC imaging
    #
    #########################################################################
    recipe = stimela.Recipe("1GC Imaging Engine", ms_dir=MSDIR)
    # imaging
    for target_field, target in zip(target_fields, targets):
        imname = cfg.obs.basename + "_1GC_" + target.name
        recipe.add("cab/wsclean", "wsclean_%d" % field_index[target.name],
            {
                "msname"            : cfg.obs.msfile,
                "column"            : cfg.wsclean_image.column,
                "weight"            : "briggs %.2f"%(cfg.wsclean_image.robust),
                "npix"              : cfg.obs.im_npix,
                "cellsize"          : cfg.obs.angular_resolution*cfg.obs.sampling,
                "clean_iterations"  : cfg.wsclean_image.clean_iterations,
                "mgain"             : cfg.wsclean_image.mgain,
                "channelsout"       : cfg.obs.im_numchans,
                "joinchannels"      : cfg.wsclean_image.joinchannels,
                "field"             : target_field,
                "name"              : imname,
                "auto-threshold"    : cfg.wsclean_image.autothreshold,
            },
            input=INPUT, output=OUTPUT,
            label="image_%d::wsclean" % target_field)

    recipe.run(["image_%d" % t for t in target_fields])

    #########################################################################
    #
    # Phase only self calibration (2nd gen)
    #
    #########################################################################
    recipe = stimela.Recipe("Initial phase-only 2GC Pipeline", ms_dir=MSDIR)

    # Add bitflag column. To keep track of flagsets
    recipe.add("cab/msutils", "msutils",
        {
            'command'    : 'prep',
            'msname'     : cfg.obs.msfile,
        },
        input=INPUT, output=OUTPUT,
        label="prepms::Adds flagsets")


    # Copy CORRECTED_DATA to DATA, so we can start selfcal
    recipe.add("cab/msutils", "shift_columns",
        {
            "command"           : cfg.move_corrdata_to_data.command,
            "msname"            : cfg.obs.msfile,
            "fromcol"           : cfg.move_corrdata_to_data.fromcol,
            "tocol"             : cfg.move_corrdata_to_data.tocol,
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
                "image"         : "%s:output"%imname_mfs,
                "outfile"       : model_prefix+".fits",
                "thresh_pix"    : cfg.source_find0.thresh_pix,
                "thresh_isl"    : cfg.source_find0.thresh_isl,
                "port2tigger"   :   True,
            },
            input=INPUT, output=OUTPUT,
            label="source_find0_%d:: Extract sources from previous round of cal" % target_field)

        # Stitch wsclean channel images into a cube
        cubename = cfg.obs.basename + "_1GC_" + plot_name[target.name] + "-CLEAN_cube.fits"
        recipe.add("cab/fitstool", "fitstool",
            {
                "image"     : [ '%s-%04d-image.fits:output'%(imname_prefix, a) for a in range(cfg.obs.im_numchans) ],
                "output"    : cubename,
                "stack"     : True,
                "fits-axis" : 3,
            },
            input=INPUT, output=OUTPUT,
            label="stitch_cube0_%d:: Stitch MFS image slices into a cube" % target_field)

        # Add SPIs
        recipe.add("cab/specfit", "add_SPIs_LSM0",
            {
                "image"                 :   "%s:output"%cubename,
                "output-spi-image"      :   "%s-spi.fits"%imname_prefix,
                "output-spi-error-image":   "%s-spi.error.fits"%imname_prefix,
                "input-skymodel"        :   "%s:output"%model_name, # model to which SPIs must be added
                "output-skymodel"       :   model_name,
                "tolerance-range"       :   cfg.specfit0.tol,
                "freq0"                 :   cfg.obs.freq_0,    # reference frequency for SPI calculation
                "sigma-level"           :   cfg.specfit0.sigma,
            },
            input=INPUT, output=OUTPUT,
            label="SPI0_%d::Add SPIs to LSM" % target_field)

        # Selfcal and subtract brightest sources
        recipe.add("cab/calibrator", "Initial_Gjones_subtract_LSM0",
            {
                "skymodel"  :   "%s:output"%model_name,
                "label"     :   cfg.selfcal0.label,
                "msname"    :   cfg.obs.msfile,
                "threads"   :   cfg.selfcal0.ncpu,
                "column"    :   cfg.selfcal0.column,
                "output-data"    :   cfg.selfcal0.output,
                "Gjones"    :   cfg.selfcal0.gjones,
                "Gjones-solution-intervals" : [int(math.ceil(float(cfg.obs.gain_sol_int[:-1]))),
                                               int(cfg.obs.nchans /
                                                   float(100))],
                "DDjones-smoothing-intervals" :  cfg.selfcal0.ddjones_smoothing,
                # TODO: MeerKAT beams need to go in this section
                "Ejones"    :   cfg.selfcal0.ejones,
                "beam-files-pattern" : cfg.selfcal0.beam_files_pattern,
                "beam-l-axis" : cfg.selfcal0.beam_l_axis,
                "beam-m-axis" : cfg.selfcal0.beam_m_axis,
                "Gjones-ampl-clipping"  :   cfg.selfcal0.gjones_ampl_clipping,
                "Gjones-ampl-clipping-low"  :   cfg.selfcal0.gjones_ampl_clipping_low,
                "Gjones-ampl-clipping-high"  :   cfg.selfcal0.gjones_ampl_clipping_high,
                "Gjones-matrix-type" : cfg.selfcal0.gjones_matrix_type,
                "make-plots" : True,
                "field-id" : target_field
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
                "image"     :   '%s:output'%imname_mfs,
                "output"    :   maskname,
                "sigma"     :   cfg.cleanmask0.sigma,
                "iters"     :   cfg.cleanmask0.iters,
                "boxes"    :   cfg.cleanmask0.kernel,
            },
            input=INPUT, output=OUTPUT,
            label="MSK_SC0_%d::Make clean mask" % target_field)

        # make another mfs image
        imname_prefix = cfg.obs.basename + "_SC0_" + plot_name[target.name]
        imname_mfs = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC0_%d" % target_field,
            {
                "msname"            : cfg.obs.msfile,
                "column"            : cfg.wsclean_selfcal0.column,
                "weight"            : "briggs %2.f"%(cfg.wsclean_selfcal0.robust),
                "npix"              : cfg.obs.im_npix,
                "cellsize"          : cfg.obs.angular_resolution*cfg.obs.sampling,
                "clean_iterations"  : cfg.wsclean_selfcal0.clean_iterations,
                "mgain"             : cfg.wsclean_selfcal0.mgain,
                "channelsout"       : cfg.obs.im_numchans,
                "joinchannels"      : cfg.wsclean_selfcal0.joinchannels,
                "field"             : str(target_field),
                "name"              : imname_prefix,
                "fitsmask"          : '%s:output' % maskname,
                "auto-threshold"    : cfg.wsclean_selfcal0.autothreshold,
            },
            input=INPUT, output=OUTPUT,
            label="image_SC0_%d::wsclean" % target_field)

        # full restore into the mfs
        imname_mfs_fullrest = imname_prefix + "-MFS-image.full_restore.fits"
        recipe.add("cab/tigger_restore", "full_restore_SC0_%d" % target_field,
            {
                "input-image"       : "%s:output" % imname_mfs,
                "input-skymodel"    : "%s:output" % model_name,
                "output-image"      : '%s:output' % imname_mfs_fullrest,
                "freq"              : cfg.obs.freq_0,
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

    #########################################################################
    #
    # Phase Amplitude self calibration (2nd gen)
    #
    # After most of the significant phase error has been corrected
    # we should be able to clean deeper, so rerun this time with amplitude
    #########################################################################
    recipe = stimela.Recipe("Phase amplitude 2GC Pipeline", ms_dir=MSDIR)

    # Initial selfcal loop
    for target_field, target in zip(target_fields, targets):
        # Extract sources in mfs clean image to build initial sky model
        imname_prefix = cfg.obs.basename + "_SC0_" + plot_name[target.name]
        imname_mfs = imname_prefix + "-MFS-image.fits"

        # Construct LSM from sky cleaned after first SC subtraction
        old_model_prefix = cfg.obs.basename + "_LSM0_" + plot_name[target.name]
        old_model_name = old_model_prefix + ".lsm.html"
        model_prefix = cfg.obs.basename + "_LSM1_" + plot_name[target.name]
        model_name = model_prefix + ".lsm.html"
        model_master = cfg.obs.basename + \
                       "_LSM_concatenated_" + \
                       plot_name[target.name] + \
                       ".lsm.html"
        recipe.add("cab/pybdsm", "extract_sources_%d" % target_field,
            {
                "image"         : "%s:output"%imname_mfs,
                "outfile"       : model_prefix+".fits",
                "thresh_pix"    : cfg.source_find1.thresh_pix,
                "thresh_isl"    : cfg.source_find1.thresh_isl,
                "port2tigger"   :   True,
            },
            input=INPUT, output=OUTPUT,
            label="source_find1_%d:: Extract sources from previous round of cal" % target_field)

        # Stitch wsclean channel images of the previous run into a cube
        cubename = cfg.obs.basename + "_SC0_" + plot_name[target.name] + "-CLEAN_cube.fits"
        recipe.add("cab/fitstool", "fitstool",
            {
                "image"     : [ '%s-%04d-image.fits:output'%(imname_prefix, a) for a in range(cfg.obs.im_numchans) ],
                "output"    : cubename,
                "stack"     : True,
                "fits-axis" : 3,
            },
            input=INPUT, output=OUTPUT,
            label="stitch_cube1_%d:: Stitch MFS image slices into a cube" % target_field)

        # Add SPIs
        recipe.add("cab/specfit", "add_SPIs_LSM0",
            {
                "image"                 :   "%s:output"%cubename,
                "output-spi-image"      :   "%s-spi.fits"%imname_prefix,
                "output-spi-error-image":   "%s-spi.error.fits"%imname_prefix,
                "input-skymodel"        :   "%s:output"%model_name, # model to which SPIs must be added
                "output-skymodel"       :   model_name,
                "tolerance-range"       :   cfg.specfit1.tol,
                "freq0"                 :   cfg.obs.freq_0,    # reference frequency for SPI calculation
                "sigma-level"           :   cfg.specfit1.sigma,
            },
            input=INPUT, output=OUTPUT,
            label="SPI1_%d::Add SPIs to LSM" % target_field)

        # Stitch LSM 0 and LSM 1 together into a master LSM file from which
        # we can do prediction for selfcal
        recipe.add("cab/tigger_convert", "stitch_LSM0_LSM1",
            {
                "input-skymodel"        :   "%s:output" % old_model_name,
                "output-skymodel"       :   "%s:output" % model_master,
                "append"                :   "%s:output" % model_name,
            },
            input=INPUT, output=OUTPUT,
            label="stitch_lsms1_%d::Create master lsm file from SC0 and SC1" %
                   target_field )

        # Selfcal and subtract brightest sources
        recipe.add("cab/calibrator", "Initial_Gjones_subtract_LSM1",
            {
                "skymodel"  :   "%s:output" % model_master,
                "label"     :   cfg.selfcal1.label,
                "msname"    :   cfg.obs.msfile,
                "threads"   :   cfg.selfcal1.ncpu,
                "column"    :   cfg.selfcal1.column,
                "output-data"    :   cfg.selfcal1.output,
                "Gjones"    :   cfg.selfcal1.gjones,
                "Gjones-solution-intervals" : [int(math.ceil(float(cfg.obs.gain_sol_int[:-1]))),
                                               int(cfg.obs.nchans /
                                                   float(100))],
                "DDjones-smoothing-intervals" :  cfg.selfcal1.ddjones_smoothing,
                # TODO: MeerKAT beams need to go in this section
                "Ejones"    :   cfg.selfcal1.ejones,
                "beam-files-pattern" : cfg.selfcal1.beam_files_pattern,
                "beam-l-axis" : cfg.selfcal1.beam_l_axis,
                "beam-m-axis" : cfg.selfcal1.beam_m_axis,
                "Gjones-ampl-clipping"  :   cfg.selfcal1.gjones_ampl_clipping,
                "Gjones-ampl-clipping-low"  :   cfg.selfcal1.gjones_ampl_clipping_low,
                "Gjones-ampl-clipping-high"  :   cfg.selfcal1.gjones_ampl_clipping_high,
                "Gjones-matrix-type" : cfg.selfcal1.gjones_matrix_type,
                "make-plots" : True,
                "field-id" : target_field
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

        # create a mask for this round of selfcal
        maskname = imname_prefix + "_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask",
            {
                "image"     :   '%s:output'%imname_mfs,
                "output"    :   maskname,
                "sigma"     :   cfg.cleanmask1.sigma,
                "iters"     :   cfg.cleanmask1.iters,
                "boxes"    :   cfg.cleanmask1.kernel,
            },
            input=INPUT, output=OUTPUT,
            label="MSK_SC1_%d::Make clean mask" % target_field)

        #make another mfs image
        imname_prefix = cfg.obs.basename + "_SC1_" + plot_name[target.name]
        recipe.add("cab/wsclean", "wsclean_SC1_%d" % target_field,
            {
                "msname"            : cfg.obs.msfile,
                "column"            : cfg.wsclean_selfcal1.column,
                "weight"            : "briggs %2.f"%(cfg.wsclean_selfcal1.robust),
                "npix"              : cfg.obs.im_npix,
                "cellsize"          : cfg.obs.angular_resolution*cfg.obs.sampling,
                "clean_iterations"  : cfg.wsclean_selfcal1.clean_iterations,
                "mgain"             : cfg.wsclean_selfcal1.mgain,
                "channelsout"       : cfg.obs.im_numchans,
                "joinchannels"      : cfg.wsclean_selfcal1.joinchannels,
                "field"             : str(target_field),
                "name"              : imname_prefix,
                "fitsmask"          : '%s:output' % maskname,
                "auto-threshold"    : cfg.wsclean_selfcal1.autothreshold,
            },
            input=INPUT, output=OUTPUT,
            label="image_SC1_%d::wsclean" % target_field)

        # full restore into the mfs
        imname_mfs_fullrest = imname_prefix + "-MFS-image.full_restore.fits"
        recipe.add("cab/tigger_restore", "full_restore_SC1_%d" % target_field,
            {
                "input-image"       : imname_mfs,
                "input-skymodel"    : model_name,
                "output-image"      : '%s:output' % imname_mfs_fullrest,
                "freq"              : cfg.obs.freq_0,
            },
            input=INPUT, output=OUTPUT,
            label="fullrestore_SC1_%d::tigger_restore" % target_field)

    for ti in targets:
        # Extract sources in mfs clean image to build initial sky model
        imname_prefix = cfg.obs.basename + "_STOKES_V_RESIDUE_" + plot_name[target.name]

        # it  is common belief that the Unverse is free of
        # stokes V for the most part. So any structure left
        # in it is probably calibration artifacts
        recipe.add("cab/wsclean", "wsclean_stokes_v",
        {
            "msname"            : cfg.obs.msfile,
            "column"            : cfg.wsclean_v_residue.column,
            "weight"            : "briggs %.2f"%(cfg.wsclean_v_residue.robust),
            "npix"              : cfg.obs.im_npix,
            "cellsize"          : cfg.obs.angular_resolution*cfg.obs.sampling,
            "clean_iterations"  : cfg.wsclean_v_residue.clean_iterations,
            "mgain"             : cfg.wsclean_v_residue.mgain,
            "channelsout"       : cfg.obs.im_numchans,
            "joinchannels"      : cfg.wsclean_v_residue.joinchannels,
            "field"             : str(target_field),
            "name"              : imname_prefix,
            "pol"               : cfg.wsclean_v_residue.pol,
        },
        input=INPUT, output=OUTPUT,
        label="image_stokesv_residue_%d::wsclean image STOKES "
              "V as diagnostic" % target_field)

    # Initial selfcal loop
    phaseamp_2gc = []
    for target_field in target_fields:
        phaseamp_2gc += ["source_find1_%d" % target_field,
                         "stitch_cube1_%d" % target_field,
                         "SPI1_%d" % target_field,
                         "stitch_lsms1_%d" % target_field,
                         "SELFCAL1_%d" % target_field,
                         "MSK_SC1_%d" % target_field,
                         "image_SC1_%d" % target_field,
                         "fullrestore_SC1_%d" % target_field,
                        ]
    # diagnostic only:
    phaseamp_2gc += ["image_stokesv_residue_%d" % t for t in target_fields]

    recipe.run(phaseamp_2gc)

