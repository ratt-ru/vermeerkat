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
import re
import glob
import datetime
import pytz
import copy
from functools import reduce
import math
import numpy as np

import stimela
import vermeerkat
import vermeerkat.caltables as vmct
import vermeerkat.config as vmc
import vermeerkat.observation as vmo
import vermeerkat.utils as vmu

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
    recipe = stimela.Recipe("Conversion Engine", ms_dir=MSDIR)

    # Convert
    recipe.add("cab/h5toms", "h5toms",
        {
            'hdf5files'  : [cfg.obs.h5file],
            'output-ms'  : cfg.obs.msfile,
            'model-data' : True,
            'full_pol'   : cfg.h5toms.full_pol,
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

    #########################################################################
    #
    # Preliminary RFI, band and autocorr flagging
    #
    #########################################################################
    recipe = stimela.Recipe("Initial flagging Engine", ms_dir=MSDIR)

    # RFI and bad channel flagging
    recipe.add("cab/rfimasker", "mask_knownrfi",
        {
            "msname" : cfg.obs.msfile,
            "mask"   : cfg.rfimask.rfi_mask_file,
        },
        input=INPUT, output=OUTPUT,
        label="rfimask::maskms")

    recipe.add("cab/autoflagger", "auto_flag_rfi",
        {
            "msname"    : cfg.obs.msfile,
            "column"    : cfg.autoflag.column,
            "strategy"  : cfg.autoflag.strategy_file,
        },
        input=INPUT, output=OUTPUT,
        label="autoflag:: Auto Flagging ms")

    recipe.add("cab/casa_flagdata", "flag_bad_start_channels",
        {
            "msname"    :   cfg.obs.msfile,
            "mode"      :   cfg.flag_bandstart.mode,
            "field"     :   cfg.flag_bandstart.field,
            "spw"       :   cfg.flag_bandstart.spw,
            "autocorr"  :   cfg.flag_bandstart.autocorr,
        },
        input=INPUT, output=OUTPUT,
        label="flag_bandstart:: Flag start of band")

    recipe.add("cab/casa_flagdata", "flag_bad_end_channels",
        {
            "msname"    :   cfg.obs.msfile,
            "mode"      :   cfg.flag_bandend.mode,
            "field"     :   cfg.flag_bandend.field,
            "spw"       :   cfg.flag_bandend.spw,
            "autocorr"  :   cfg.flag_bandend.autocorr,
        },
        input=INPUT, output=OUTPUT,
        label="flag_bandend:: Flag end of band")

    recipe.add("cab/casa_flagdata", "flag_autocorrs",
        {
            "msname"    :   cfg.obs.msfile,
            "mode"      :   cfg.flag_autocorrs.mode,
            "field"     :   cfg.flag_autocorrs.field,
            "spw"       :   cfg.flag_autocorrs.spw,
            "autocorr"  :   cfg.flag_autocorrs.autocorr,
        },
        input=INPUT, output=OUTPUT,
        label="flag_autocorrs:: Flag auto correlations")

    # Run flagging run!
    init_flagging = ["rfimask",
                     "autoflag",
                     "flag_bandstart",
                     "flag_bandend",
                     "flag_autocorrs"]
    recipe.run(init_flagging)

    #########################################################################
    #
    # Initial 1GC
    #
    # We first calibrate the calibrators for flagging
    # purposes, we will want to discard these solutions when we
    # did mitigation flagging and generate, hopefully, improved 1GC
    # solutions
    #########################################################################

    # Selects the flux scale based on selected bandpass calibrator
    # if it is in our southern standard then grab the coefficients
    # to plug into CASA, otherwise fall back to the CASA standard
    # as last resort and then fall over if it ain't in there either
    if bandpass_cal.name in calibrator_db:
        aghz = calibrator_db[bandpass_cal.name]["a_ghz"]
        bghz = calibrator_db[bandpass_cal.name]["b_ghz"]
        cghz = calibrator_db[bandpass_cal.name]["c_ghz"]
        dghz = calibrator_db[bandpass_cal.name]["d_ghz"]

        # Find the brightness at reference frequency
        I, a, b, c, d = vmct.convert_pb_to_casaspi(
            cfg.obs.freq_0 / 1e9 - (cfg.obs.nchans // 2) * cfg.obs.chan_bandwidth / 1e9,
            cfg.obs.freq_0 / 1e9 + (cfg.obs.nchans // 2) * cfg.obs.chan_bandwidth / 1e9,
            cfg.obs.freq_0 / 1e9, aghz, bghz, cghz, dghz)

        vermeerkat.log.info("Using bandpass calibrator %s "
                            "with brightness of %.4f Jy "
                            "(spix = [%.6f, %.6f, %.6f, %.6f]) "
                            "(@ %.2f MHz) "
                            "as the flux scale reference" %
                            (bandpass_cal.name,
                             I, a, b, c, d,
                             cfg.obs.freq_0 / 1e6))
        # 1GC Calibration
        setjy_options = {
                "msname"        :   cfg.obs.msfile,
                "field"         :   str(bpcal_field),
                "standard"      :   cfg.setjy_manual.standard,
                "fluxdensity"   :   I,
                "spix"          :   [a, b, c, d],
                "reffreq"       :   "%.2fGHz" % (cfg.obs.freq_0 / 1e9),
                "usescratch"    :   cfg.setjy_manual.usescratch,
                "scalebychan"   :   cfg.setjy_manual.scalebychan,
                "spw"           :   cfg.setjy_manual.spw,
        }
    else:
        vermeerkat.log.warn("Looks like your flux reference '%s' is not "
            "in our standard. We will now try to fall back to the CASA "
            "standard you specified in your config file. "
            "Try pulling the latest "
            "VermeerKAT or if you have done "
            "so report this issue" % bandpass_cal.name)

        # If model is not in @Ben's southern calibrators, then use CASA model.
        # If this is the case, the standard must be specified in the config file
        setjy_options = {
            "msname"    :   cfg.obs.msfile,
            "standard"  :   cfg.setjy_auto.standard,
            "field"     :   str(bpcal_field)
        }

    # @mauch points out some antenna positions may
    # be slightly off. This will cause a noticible
    # decorrelation on the longest baselines, so
    # better calibrate for this one (we can use the
    # bright bandpass and gain calibrator sources).
    delaycal_opts = {
            "msname"        : cfg.obs.msfile,
            "field"         : str(gaincal_field),
            "gaintype"      : cfg.delaycal.gaintype,
            "solint"        : cfg.delaycal.solint,
            "minsnr"        : cfg.delaycal.minsnr,
            "refant"        : cfg.obs.refant,
            "caltable"      : cfg.obs.delaycal_table,
            "calmode"       : cfg.delaycal.calmode,
    }

    # The bandpass calibrator may vary in phase over time
    # so lets do a preliminary phase cal to correct for this
    bp_phasecal_opts = {
            "msname"        :   cfg.obs.msfile,
            "caltable"      :   cfg.obs.phasecal_table,
            "field"         :   str(bpcal_field),
            "refant"        :   cfg.obs.refant,
            "calmode"       :   cfg.phase0.calmode,
            "solint"        :   cfg.obs.gain_sol_int,
            "solnorm"       :   cfg.phase0.solnorm,
            "minsnr"        :   cfg.phase0.minsnr,
            "gaintable"     :   [cfg.obs.delaycal_table],
    }

    # Then a bandpass calibration, lets work out a solution
    # per scan interval to see if things vary signifcantly
    # over time... this will be very bad for your reduction.
    bandpass_cal_opts =  {
            "msname"        :   cfg.obs.msfile,
            "caltable"      :   cfg.obs.bpasscal_table,
            "field"         :   str(bpcal_field),
            "spw"           :   cfg.bandpass.spw,
            "refant"        :   cfg.obs.refant,
            "solnorm"       :   cfg.bandpass.solnorm,
            "combine"       :   cfg.bandpass.combine,
            "solint"        :   str(bpcal_sol_int) + "s",
            "bandtype"      :   cfg.bandpass.bandtype,
            "minblperant"   :   cfg.bandpass.minblperant,
            "minsnr"        :   cfg.bandpass.minsnr,
            "gaintable"     :   [cfg.obs.delaycal_table,
                                 cfg.obs.phasecal_table],
            "interp"        :   cfg.bandpass.interp,
    }

    # Finally we do a second order correction on the gain
    # cal source that is closest to the target fields
    gain_cal_opts = {
            "msname"       :   cfg.obs.msfile,
            "caltable"     :   cfg.obs.ampcal_table,
            "field"        :   ",".join([str(x) for
                                        x in [bpcal_field, gaincal_field]]),
            "spw"          :   cfg.gaincal.spw,
            "solint"       :   cfg.obs.gain_sol_int,
            "refant"       :   cfg.obs.refant,
            "gaintype"     :   cfg.gaincal.gaintype,
            "calmode"      :   cfg.gaincal.calmode,
            "solnorm"      :   cfg.gaincal.solnorm,
            "gaintable"    :   [cfg.obs.delaycal_table,
                                cfg.obs.phasecal_table,
                                cfg.obs.bpasscal_table],
            "interp"       :   cfg.gaincal.interp,
    }

    # Scale the gaincal solutions amplitude to that of the
    # bandpass calibrator (the flux scale reference)
    flux_scale_opts = {
            "msname"        :   cfg.obs.msfile,
            "caltable"      :   cfg.obs.ampcal_table,
            "fluxtable"     :   cfg.obs.fluxcal_table,
            "reference"     :   str(bpcal_field),
            "transfer"      :   str(gaincal_field),
            "incremental"   :   cfg.fluxscale.incremental,
    }

    # Apply gain solutions to all fields including
    # the calibrators so that we can diagnose problems more
    # easily. Later steps depend on this so don't remove
    # the gain application to calibrators
    apply_cal_opts = {
            "msname"        :   cfg.obs.msfile,
            "field"         :   ",".join([str(x) for x in
                                    target_fields + [bpcal_field, gaincal_field]]),
            "gaintable"     :   [cfg.obs.delaycal_table,
                                 cfg.obs.phasecal_table,
                                 cfg.obs.bpasscal_table,
                                 cfg.obs.fluxcal_table],
            "gainfield"     :   [str(x) for x in
                                    gaincal_field,
                                    bpcal_field,
                                    bpcal_field,
                                    gaincal_field],
            "interp"        :   cfg.applycal.interp,
            "spwmap"        :   cfg.applycal.spwmap,
            "parang"        :   cfg.applycal.parang,
    }

    recipe = stimela.Recipe("1.1GC Engine", ms_dir=MSDIR)

    recipe.add("cab/casa_setjy", "init_flux_scaling",
               setjy_options,
               input=INPUT, output=OUTPUT,
               label="setjy:: Initial flux density scaling")

    recipe.add("cab/casa_gaincal", "delay_cal",
        delaycal_opts,
        input=INPUT, output=OUTPUT,
        label="delaycal:: Delay calibration")

    recipe.add("cab/casa_gaincal", "init_phase_cal",
        bp_phasecal_opts,
        input=INPUT, output=OUTPUT,
        label="phase0:: Initial phase calibration")

    recipe.add("cab/casa_bandpass", "bandpass_cal",
        bandpass_cal_opts,
        input=INPUT, output=OUTPUT,
        label="bandpass:: Bandpass calibration")

    recipe.add("cab/casa_gaincal", "main_gain_calibration",
        gain_cal_opts,
        input=INPUT, output=OUTPUT,
        label="gaincal:: Gain calibration")

    recipe.add("cab/casa_fluxscale", "casa_fluxscale",
        flux_scale_opts,
        input=INPUT, output=OUTPUT,
        label="fluxscale:: Setting Fluxscale")

    recipe.add("cab/casa_applycal", "apply_calibration",
        apply_cal_opts,
        input=INPUT, output=OUTPUT,
        label="applycal:: Apply calibration solutions to target")

    #Go initial 1GC and further flagging
    init_1gc = [ "setjy",
                 "delaycal",
                 "phase0",
                 "bandpass",
                 "gaincal",
                 "fluxscale",
                 "applycal"]
    recipe.run(init_1gc)

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
    recipe.add("cab/politsiyakat", "flag_malfunctioning_antennas",
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
        label="flag_baseline_phases_bp:: Flag baselines based on calibrator phases")

    post_1gc_flagging = [ "autoflag_corrected_vis",
                          "flag_baseline_phases_bp"]
    recipe.run(post_1gc_flagging)

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
    recipe = stimela.Recipe("1.2GC Engine", ms_dir=MSDIR)

    recipe.add("cab/casa_setjy", "init_flux_scaling_postmit",
               setjy_options,
               input=INPUT, output=OUTPUT,
               label="setjy:: Initial flux density scaling")

    delaycal_opts["caltable"] = cfg.obs.post_mit_delaycal_table
    recipe.add("cab/casa_gaincal", "delay_cal",
        delaycal_opts,
        input=INPUT, output=OUTPUT,
        label="delaycal:: Delay calibration")

    bp_phasecal_opts["caltable"] = cfg.obs.post_mit_phasecal_table
    bp_phasecal_opts["gaintable"] = [cfg.obs.post_mit_delaycal_table]
    recipe.add("cab/casa_gaincal", "init_phase_cal_postmit",
        bp_phasecal_opts,
        input=INPUT, output=OUTPUT,
        label="phase0:: Initial phase calibration")

    bandpass_cal_opts["caltable"] = cfg.obs.post_mit_bpasscal_table
    bandpass_cal_opts["gaintable"] = [cfg.obs.post_mit_delaycal_table,
                                cfg.obs.post_mit_phasecal_table]
    recipe.add("cab/casa_bandpass", "bandpass_cal_postmit",
        bandpass_cal_opts,
        input=INPUT, output=OUTPUT,
        label="bandpass:: Bandpass calibration")

    gain_cal_opts["caltable"] = cfg.obs.post_mit_ampcal_table
    gain_cal_opts["gaintable"] = [cfg.obs.post_mit_delaycal_table,
                            cfg.obs.post_mit_phasecal_table,
                            cfg.obs.post_mit_bpasscal_table]
    recipe.add("cab/casa_gaincal", "main_gain_calibration_postmit",
        gain_cal_opts,
        input=INPUT, output=OUTPUT,
        label="gaincal:: Gain calibration")

    flux_scale_opts["caltable"] = cfg.obs.post_mit_ampcal_table
    flux_scale_opts["fluxtable"] = cfg.obs.post_mit_fluxcal_table
    recipe.add("cab/casa_fluxscale", "casa_fluxscale_postmit",
        flux_scale_opts,
        input=INPUT, output=OUTPUT,
        label="fluxscale:: Setting Fluxscale")

    apply_cal_opts["gaintable"] =[cfg.obs.post_mit_delaycal_table,
                             cfg.obs.post_mit_phasecal_table,
                             cfg.obs.post_mit_bpasscal_table,
                             cfg.obs.post_mit_fluxcal_table]
    recipe.add("cab/casa_applycal", "apply_calibration_postmit",
        apply_cal_opts,
        input=INPUT, output=OUTPUT,
        label="applycal:: Apply calibration solutions to target")

    second_1gc = [ "setjy",
                   "delaycal",
                   "phase0",
                   "bandpass",
                   "gaincal",
                   "fluxscale",
                   "applycal",
                 ]
    try:
        recipe.run(second_1gc)
        pass
    except:
        vermeerkat.log.warning("Post mitigation 1GC calibration failure. "
                               "Inspect the logs. We will continue with "
                               "the solutions derived pre-mitigation "
                               "flagging, but with mitigation flags applied")

    #########################################################################
    #
    # Post 1GC diagnostic plots
    #
    # With our best solutions in hand we plot some results for the observer
    #########################################################################
    recipe = stimela.Recipe("1GC Diagnostics Engine", ms_dir=MSDIR)

    # Diagnostic: amplitude vs uv dist of the bandpass calibrator
    # @mauch points out we expect all baselines to observe the same amplitude
    # for the point source-like bandpass calibrator
    plot_ampuvdist_opts = {
            "msname"            : cfg.obs.msfile,
            "xaxis"             : cfg.plot_ampuvdist.xaxis,
            "yaxis"             : cfg.plot_ampuvdist.yaxis,
            "xdatacolumn"       : cfg.plot_ampuvdist.xdatacolumn,
            "ydatacolumn"       : cfg.plot_ampuvdist.ydatacolumn,
            "field"             : str(bpcal_field),
            "iteraxis"          : cfg.plot_ampuvdist.iteraxis,
            "avgchannel"        : cfg.plot_ampuvdist.avgchannel,
            "avgtime"           : cfg.plot_ampuvdist.avgtime,
            "coloraxis"         : cfg.plot_ampuvdist.coloraxis,
            "expformat"         : cfg.plot_ampuvdist.expformat,
            "exprange"          : cfg.plot_ampuvdist.exprange,
            "plotfile"          : cfg.obs.basename + "_bp_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "ampuvdist.png",
            "overwrite"         :   True,
    }

    #Diagnostic: phase vs uv dist of the bandpass calibrator
    plot_phaseuvdist_opts = {
            "msname"            : cfg.obs.msfile,
            "xaxis"             : cfg.plot_phaseuvdist.xaxis,
            "yaxis"             : cfg.plot_phaseuvdist.yaxis,
            "xdatacolumn"       : cfg.plot_phaseuvdist.xdatacolumn,
            "ydatacolumn"       : cfg.plot_phaseuvdist.ydatacolumn,
            "field"             : str(bpcal_field),
            "iteraxis"          : cfg.plot_phaseuvdist.iteraxis,
            "avgchannel"        : cfg.plot_phaseuvdist.avgchannel,
            "avgtime"           : cfg.plot_phaseuvdist.avgtime,
            "coloraxis"         : cfg.plot_phaseuvdist.coloraxis,
            "expformat"         : cfg.plot_phaseuvdist.expformat,
            "exprange"          : cfg.plot_phaseuvdist.exprange,
            "plotfile"          : cfg.obs.basename + "_bp_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "phaseuvdist.png",
            "overwrite"         :   True,
    }

    # Diagnostic: amplitude vs phase of bp calibrator per antenna
    plot_phaseball_opts = {
            "msname"            : cfg.obs.msfile,
            "xaxis"             : cfg.plot_phaseball.xaxis,
            "yaxis"             : cfg.plot_phaseball.yaxis,
            "xdatacolumn"       : cfg.plot_phaseball.xdatacolumn,
            "ydatacolumn"       : cfg.plot_phaseball.ydatacolumn,
            "field"             : str(bpcal_field),
            "iteraxis"          : cfg.plot_phaseball.iteraxis,
            "avgchannel"        : cfg.plot_phaseball.avgchannel,
            "avgtime"           : cfg.plot_phaseball.avgtime,
            "coloraxis"         : cfg.plot_phaseball.coloraxis,
            "expformat"         : cfg.plot_phaseball.expformat,
            "exprange"          : cfg.plot_phaseball.exprange,
            "plotfile"          : cfg.obs.basename + "_bp_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "phaseball.png",
            "overwrite"         :   True,
    }

    # Diagnostic: amplitude vs frequency of bp calibrator
    plot_amp_freq_opts = {
            "msname"            : cfg.obs.msfile,
            "xaxis"             : cfg.plot_amp_freq.xaxis,
            "yaxis"             : cfg.plot_amp_freq.yaxis,
            "xdatacolumn"       : cfg.plot_amp_freq.xdatacolumn,
            "ydatacolumn"       : cfg.plot_amp_freq.ydatacolumn,
            "field"             : str(bpcal_field),
            "iteraxis"          : cfg.plot_amp_freq.iteraxis,
            "avgchannel"        : cfg.plot_amp_freq.avgchannel,
            "avgtime"           : cfg.plot_amp_freq.avgtime,
            "coloraxis"         : cfg.plot_amp_freq.coloraxis,
            "expformat"         : cfg.plot_amp_freq.expformat,
            "exprange"          : cfg.plot_amp_freq.exprange,
            "plotfile"          : cfg.obs.basename + "_bp_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "band.png",
            "overwrite"         :   True,
    }

    # Diagnostic: phase vs time of bp calibrator
    # @oms points out slopes in this will indicate problems with
    # digitizer reference timing / probably also any uncorrected
    # antenna positions
    plot_phase_time_opts = {
            "msname"            : cfg.obs.msfile,
            "xaxis"             : cfg.plot_phase_time.xaxis,
            "yaxis"             : cfg.plot_phase_time.yaxis,
            "xdatacolumn"       : cfg.plot_phase_time.xdatacolumn,
            "ydatacolumn"       : cfg.plot_phase_time.ydatacolumn,
            "field"             : str(bpcal_field),
            "iteraxis"          : cfg.plot_phase_time.iteraxis,
            "avgchannel"        : cfg.plot_phase_time.avgchannel,
            "avgtime"           : cfg.plot_phase_time.avgtime,
            "coloraxis"         : cfg.plot_phase_time.coloraxis,
            "expformat"         : cfg.plot_phase_time.expformat,
            "exprange"          : cfg.plot_phase_time.exprange,
            "plotfile"          : cfg.obs.basename + "_bp_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "phasevtime.png",
            "overwrite"         :   True,
    }
    # Diagnostic: phase vs time of bp calibrator
    # For similar purposes as phase vs freq
    plot_phase_freq_opts = {
            "msname"            : cfg.obs.msfile,
            "xaxis"             : cfg.plot_phase_freq.xaxis,
            "yaxis"             : cfg.plot_phase_freq.yaxis,
            "xdatacolumn"       : cfg.plot_phase_freq.xdatacolumn,
            "ydatacolumn"       : cfg.plot_phase_freq.ydatacolumn,
            "field"             : str(bpcal_field),
            "iteraxis"          : cfg.plot_phase_freq.iteraxis,
            "avgchannel"        : cfg.plot_phase_freq.avgchannel,
            "avgtime"           : cfg.plot_phase_freq.avgtime,
            "coloraxis"         : cfg.plot_phase_freq.coloraxis,
            "expformat"         : cfg.plot_phase_freq.expformat,
            "exprange"          : cfg.plot_phase_freq.exprange,
            "plotfile"          : cfg.obs.basename + "_bp_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "phasevfreq.png",
            "overwrite"         :   True,
    }

    recipe.add("cab/casa_plotms", "plot_amp_v_uv_dist_of_bp_calibrator",
        copy.deepcopy(plot_ampuvdist_opts),
        input=INPUT, output=OUTPUT,
        label="plot_ampuvdist_bp:: Diagnostic plot of amplitude with uvdist")

    recipe.add("cab/casa_plotms", "plot_phase_v_uv_dist_of_bp_calibrator",
        copy.deepcopy(plot_phaseuvdist_opts),
        input=INPUT, output=OUTPUT,
        label="plot_phaseuvdist_bp:: Diagnostic plot of phase with uvdist")

    recipe.add("cab/casa_plotms", "plot_amp_v_phase_of_bp_calibrator",
        copy.deepcopy(plot_phaseball_opts),
        input=INPUT, output=OUTPUT,
        label="plot_phaseball_bp:: Diagnostic plot of phaseball")

    recipe.add("cab/casa_plotms", "plot_amp_v_freq_of_bp_calibrator",
        copy.deepcopy(plot_amp_freq_opts),
        input=INPUT, output=OUTPUT,
        label="plot_amp_freq_bp:: Diagnostic plot of band")

    recipe.add("cab/casa_plotms", "plot_phase_vs_time_of_bp_calibrator",
        copy.deepcopy(plot_phase_time_opts),
        input=INPUT, output=OUTPUT,
        label="plot_phase_time_bp:: Diagnostic plot of phase with time")

    recipe.add("cab/casa_plotms", "plot_phase_vs_freq_of_bp_calibrator",
        copy.deepcopy(plot_phase_freq_opts),
        input=INPUT, output=OUTPUT,
        label="plot_phase_freq_bp:: Diagnostic plot of phase with freq")

    plot_ampuvdist_opts["field"] = str(gaincal_field)
    plot_ampuvdist_opts["plotfile"] = cfg.obs.basename + "_gc_" +\
                                      plot_name[gain_cal.name] + "_" +\
                                      "ampuvdist.png"
    recipe.add("cab/casa_plotms", "plot_amp_v_uv_dist_of_gc_calibrator",
        plot_ampuvdist_opts,
        input=INPUT, output=OUTPUT,
        label="plot_ampuvdist_gc:: Diagnostic plot of amplitude with uvdist")

    plot_phaseuvdist_opts["field"] = str(gaincal_field)
    plot_phaseuvdist_opts["plotfile"] = cfg.obs.basename + "_gc_" +\
                                      plot_name[gain_cal.name] + "_" +\
                                      "phaseuvdist.png"
    recipe.add("cab/casa_plotms", "plot_phase_v_uv_dist_of_gc_calibrator",
        plot_phaseuvdist_opts,
        input=INPUT, output=OUTPUT,
        label="plot_phaseuvdist_gc:: Diagnostic plot of phase with uvdist")

    plot_phaseball_opts["field"] = str(gaincal_field)
    plot_phaseball_opts["plotfile"] = cfg.obs.basename + "_gc_" +\
                                      plot_name[gain_cal.name] + "_" +\
                                      "phaseball.png"
    recipe.add("cab/casa_plotms", "plot_amp_v_phase_of_gc_calibrator",
        plot_phaseball_opts,
        input=INPUT, output=OUTPUT,
        label="plot_phaseball_gc:: Diagnostic plot of phaseball")

    plot_amp_freq_opts["field"] = str(gaincal_field)
    plot_amp_freq_opts["plotfile"] = cfg.obs.basename + "_gc_" +\
                                      plot_name[gain_cal.name] + "_" +\
                                      "band.png"
    recipe.add("cab/casa_plotms", "plot_amp_v_freq_of_gc_calibrator",
        plot_amp_freq_opts,
        input=INPUT, output=OUTPUT,
        label="plot_amp_freq_gc:: Diagnostic plot of band")

    plot_phase_time_opts["field"] = str(gaincal_field)
    plot_phase_time_opts["plotfile"] = cfg.obs.basename + "_gc_" +\
                                       plot_name[gain_cal.name] + "_" +\
                                       "phasevtime.png"

    plot_phase_time_opts["field"] = str(gaincal_field)
    plot_phase_time_opts["plotfile"] = cfg.obs.basename + "_gc_" +\
                                      plot_name[gain_cal.name] + "_" +\
                                      "phasevtime.png"
    recipe.add("cab/casa_plotms", "plot_phase_vs_time_of_gc_calibrator",
        plot_phase_time_opts,
        input=INPUT, output=OUTPUT,
        label="plot_phase_time_gc:: Diagnostic plot of phase with time")

    plot_phase_freq_opts["field"] = str(gaincal_field)
    plot_phase_freq_opts["plotfile"] = cfg.obs.basename + "_gc_" +\
                                      plot_name[gain_cal.name] + "_" +\
                                      "phasevfreq.png"
    recipe.add("cab/casa_plotms", "plot_phase_vs_freq_of_gc_calibrator",
        plot_phase_freq_opts,
        input=INPUT, output=OUTPUT,
        label="plot_phase_freq_gc:: Diagnostic plot of phase with freq")

    # Diagnostic only: image bandpass
    recipe.add("cab/wsclean", "wsclean_bandpass",
        {
            "msname"            : cfg.obs.msfile,
            "column"            : cfg.wsclean_bandpass.column,
            "weight"            : "briggs %.2f"%(cfg.wsclean_bandpass.robust),
            "npix"              : cfg.obs.im_npix,
            "cellsize"          : cfg.obs.angular_resolution * cfg.obs.sampling,
            "clean_iterations"  : cfg.wsclean_bandpass.clean_iterations,
            "mgain"             : cfg.wsclean_bandpass.mgain,
            "channelsout"       : cfg.obs.im_numchans,
            "joinchannels"      : cfg.wsclean_bandpass.joinchannels,
            "field"             : str(bpcal_field),
            "name"              : cfg.obs.basename + "_bp_" + plot_name[bandpass_cal.name],
        },
        input=INPUT, output=OUTPUT,
        label="image_bandpass::wsclean")

    # Diagnostic only: image gaincal
    recipe.add("cab/wsclean", "wsclean_gain",
        {
            "msname"            : cfg.obs.msfile,
            "column"            : cfg.wsclean_gain.column,
            "weight"            : "briggs %.2f"%(cfg.wsclean_gain.robust),
            "npix"              : cfg.obs.im_npix,
            "cellsize"          : cfg.obs.angular_resolution*cfg.obs.sampling,
            "clean_iterations"  : cfg.wsclean_gain.clean_iterations,
            "mgain"             : cfg.wsclean_gain.mgain,
            "channelsout"       : cfg.obs.im_numchans,
            "joinchannels"      : cfg.wsclean_gain.joinchannels,
            "field"             : str(gaincal_field),
            "name"              : cfg.obs.basename + "_gc_" + plot_name[gain_cal.name],
        },
        input=INPUT, output=OUTPUT,
        label="image_gain::wsclean")

    diagnostics_1gc = [ "plot_ampuvdist_bp",
                        "plot_phaseuvdist_bp",
                        "plot_phaseball_bp",
                        "plot_amp_freq_bp",
                        "plot_phase_time_bp",
                        "plot_phase_freq_bp",
                        "plot_ampuvdist_gc",
                        "plot_phaseuvdist_gc",
                        "plot_phaseball_gc",
                        "plot_amp_freq_gc",
                        "plot_phase_time_gc",
                        "plot_phase_freq_gc",
                      ]
    diagnostics_1gc += ["image_bandpass", "image_gain"] # diagnostic only
    recipe.run(diagnostics_1gc)

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

