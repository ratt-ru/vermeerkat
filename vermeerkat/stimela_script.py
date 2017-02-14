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
from functools import reduce

import numpy as np

import stimela
import vermeerkat
import vermeerkat.caltables as vmct
import vermeerkat.config as vmc
import vermeerkat.observation as vmo
import vermeerkat.utils as vmu

# So that we can access GLOBALS pass through to the run command
stimela.register_globals()

# args is a global string variable representing a python list of strings
# Evaluate it to get the list of strings back
args = ast.literal_eval(args)

# Load in the configuration
cfg = vmc.configuration(args)
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

    # Categories the fields observed in each scan
    field_index, bpcals, gaincals, targets = vmu.categorise_fields(scans)

    # Alias a long name
    default_bpcal = cfg.general.bandpass_calibrator

    if default_bpcal:
        # If a default bandpass calibrator has been specified,
        # filter out other bandpass calibrators. We still need
        # to extract it from the standard below
        vermeerkat.log.info("Bandpass calibrator manually "
                            "set to '%s'. Other bandpass "
                            "calibrators will be ignored." % default_bpcal)
        bpcals = [b for b in bpcals if b.name == default_bpcal]

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
    gaincal_index, gain_cal = vmu.select_gain_calibrator(targets, gaincals)
    # Compute the observation time spent on gain calibrators
    gaincal_scan_times = vmu.total_scan_times(field_scan_map, gaincals)

    # Remove gain calibrator from the bandpass calibrators
    bpcals = [x for x in bpcals if not x.name == gain_cal.name]

    # Compute observation time on bandpass calibrators
    bandpass_scan_times = vmu.total_scan_times(field_scan_map, bpcals)
    # Select the bandpass calibrator
    bpcal_index, bandpass_cal = vmu.select_bandpass_calibrator(bpcals,
                                                bandpass_scan_times)

    # Choose the solution interval for the bandpass calibrator
    # by choosing the bandpass calibrator's minimum scan length
    bpcal_sol_int = min(s.length for s in field_scan_map[bandpass_cal.name])

    # Compute observation time on target
    target_scan_times = vmu.total_scan_times(field_scan_map, targets)

    # Get field ids for bandpass, gaincal and target
    bpcal_field = field_index[bandpass_cal.name]
    gaincal_field = field_index[gain_cal.name]
    target_fields = [field_index[t.name] for t in targets]

    # Log useful information
    vermeerkat.log.info("The following fields were observed:")

    # Sort fields by index
    for k, v in sorted(field_index.items(), key=lambda (k, v): v):
        tags = set.union(*(set(s.tags) for s in field_scan_map[k]))
        scan_seconds = sum(s.length for s in field_scan_map[k])
        vermeerkat.log.info("\t %d: %s %s %s" % (v, k.ljust(20),
                    vmu.fmt_seconds(scan_seconds).ljust(20),
                    [t for t in tags]))

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
                             vmu.fmt_seconds(scan_seconds)))
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

    # All good to go fire up the pipeline
    recipe = stimela.Recipe("1GC Pipeline", ms_dir=MSDIR)

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


    # RFI and bad channel flagging
    recipe.add("cab/rfimasker", "mask_stuff",
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

    # Custom user specified flags
    # Casa does not sensibly do nothing if the defaults are specified
    # removing this step
#    userflags = {"msname" : cfg.obs.msfile}
#    if cfg.flag_userflags.mode is not None:
#        userflags["mode"] = cfg.flag_userflags.mode
#    if cfg.flag_userflags.field is not None:
#        userflags["field"] = cfg.flag_userflags.field
#    if cfg.flag_userflags.antenna is not None:
#        userflags["antenna"] = cfg.flag_userflags.antenna
#    if cfg.flag_userflags.scan is not None:
#        userflags["scan"] = cfg.flag_userflags.scan
#    if cfg.flag_userflags.timerange is not None:
#        userflags["timerange"] = cfg.flag_userflags.timerange
#    if cfg.flag_userflags.spw is not None:
#        userflags["spw"] = cfg.flag_userflags.spw
#    if cfg.flag_userflags.autocorr is not None:
#        userflags["autocorr"] = cfg.flag_userflags.autocorr
#
#    recipe.add("cab/casa_flagdata", "flag_userflags",
#        userflags,
#        input=INPUT, output=OUTPUT,
#        label="flag_userflags:: Specify any additional flags from user")

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

    #Run SetJY with our database of southern calibrators
    if bandpass_cal.name in calibrator_db:
        vermeerkat.log.warn("Looks like your flux reference '%s' is not "
                           "in our standard. Try pulling the latest "
                           "VermeerKAT or if you have done "
                           "so report this issue" % bandpass_cal.name)

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
        recipe.add("cab/casa_setjy", "init_flux_scaling",
            {
                "msname"        :   cfg.obs.msfile,
                "field"         :   str(bpcal_field),
                "standard"      :   cfg.setjy_manual.standard,
                "fluxdensity"   :   I,
                "spix"          :   [a, b, c, d],
                "reffreq"       :   "%.2fGHz" % (cfg.obs.freq_0 / 1e9),
                "usescratch"    :   cfg.setjy_manual.usescratch,
                "scalebychan"   :   cfg.setjy_manual.scalebychan,
                "spw"           :   cfg.setjy_manual.spw,
            },
            input=INPUT, output=OUTPUT,
            label="setjy:: Initial flux density scaling")
    else:
        # If model is not in @Ben's southern calibrators, then use CASA model.
        # If this is the case, the standard must be specified in the config file
        recipe.add('cab/casa_setjy', 'flux_scaling', {
            "msname"    :   cfg.obs.msfile,
            "standard"  :   cfg.setjy_auto.standard,
            "field"     :   str(bpcal_field)
            },
            input=INPUT, output=OUTPUT,
            label="setjy:: Initial flux density scaling")

    # @mauch points out some antenna positions may
    # be slightly off. This will cause a noticible
    # decorrelation on the longest baselines, so
    # better calibrate for this one (we can use the
    # bright bandpass and gain calibrator sources).
    recipe.add("cab/casa_gaincal", "delay_cal",
        {
            "msname"        : cfg.obs.msfile,
            "field"         : str(gaincal_field),
            "gaintype"      : cfg.delaycal.gaintype,
            "solint"        : cfg.delaycal.solint,
            "minsnr"        : cfg.delaycal.minsnr,
            "refant"        : cfg.obs.refant,
            "caltable"      : cfg.obs.delaycal_table,
            "calmode"       : cfg.delaycal.calmode,
        },
        input=INPUT, output=OUTPUT,
        label="delaycal:: Delay calibration")

    # The bandpass calibrator may vary in phase over time
    # so lets do a preliminary phase cal to correct for this
    recipe.add("cab/casa_gaincal", "init_phase_cal",
        {
            "msname"        :   cfg.obs.msfile,
            "caltable"      :   cfg.obs.phasecal_table,
            "field"         :   str(bpcal_field),
            "refant"        :   cfg.obs.refant,
            "calmode"       :   cfg.phase0.calmode,
            "solint"        :   cfg.obs.gain_sol_int,
            "solnorm"       :   cfg.phase0.solnorm,
            "minsnr"        :   cfg.phase0.minsnr,
            "gaintable"     :   [cfg.obs.delaycal_table],
        },
        input=INPUT, output=OUTPUT,
        label="phase0:: Initial phase calibration")

    # Then a bandpass calibration, lets work out a solution
    # per scan interval to see if things vary signifcantly
    # over time... this will be very bad for your reduction.
    recipe.add("cab/casa_bandpass", "bandpass_cal",
        {
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
        },
        input=INPUT, output=OUTPUT,
        label="bandpass:: Bandpass calibration")

    # Finally we do a second order correction on the gain
    # cal source that is closest to the target fields
    recipe.add("cab/casa_gaincal", "main_gain_calibration",
        {
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
        },
        input=INPUT, output=OUTPUT,
        label="gaincal:: Gain calibration")

    # Scale the gaincal solutions amplitude to that of the
    # bandpass calibrator (the flux scale reference)
    recipe.add("cab/casa_fluxscale", "casa_fluxscale",
        {
            "msname"        :   cfg.obs.msfile,
            "caltable"      :   cfg.obs.ampcal_table,
            "fluxtable"     :   cfg.obs.fluxcal_table,
            "reference"     :   [bandpass_cal.name],
            "transfer"      :   [gain_cal.name],
            "incremental"   :   cfg.fluxscale.incremental,
        },
        input=INPUT, output=OUTPUT,
        label="fluxscale:: Setting Fluxscale")

    # Apply gain solutions to all fields including
    # the calibrators so that we can diagnose problems more
    # easily. Later steps depend on this so don't remove
    # the gain application to calibrators
    recipe.add("cab/casa_applycal", "apply_calibration",
        {
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
        },
        input=INPUT, output=OUTPUT,
        label="applycal:: Apply calibration solutions to target")

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
    # bandpass calibrator and flag out baselines (per channel) which are not up
    # to spec

    recipe.add("cab/politsiyakat", "flag_malfunctioning_antennas",
        {
            "task"                   : cfg.flag_baseline_phases.task,
            "msname"                 : cfg.obs.msfile,
            "data_column"            : cfg.flag_baseline_phases.data_column,
            "cal_field"              : str(bpcal_field),
            "valid_phase_range"      : cfg.flag_baseline_phases.valid_phase_range,
            "max_invalid_datapoints" : cfg.flag_baseline_phases.max_invalid_datapoints,
            "output_dir"             : "",
            "nrows_chunk"            : cfg.flag_baseline_phases.nrows_chunk,
            "simulate"               : cfg.flag_baseline_phases.simulate,
        },
        input=INPUT, output=OUTPUT,
        label="flag_baseline_phases:: Flag baselines based on calibrator phases")

    # Diagnostic: amplitude vs uv dist of the bandpass calibrator
    # @mauch points out we expect all baselines to observe the same amplitude
    # for the point source-like bandpass calibrator
    recipe.add("cab/casa_plotms", "plot_amp_v_uv_dist_of_bp_calibrator",
        {
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
            "plotfile"          : cfg.obs.basename + "_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "ampuvdist.png",
            "overwrite"         :   True,
        },
        input=INPUT, output=OUTPUT,
        label="plot_ampuvdist:: Diagnostic plot of amplitude with uvdist")

    #Diagnostic: phase vs uv dist of the bandpass calibrator
    recipe.add("cab/casa_plotms", "plot_phase_v_uv_dist_of_bp_calibrator",
        {
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
            "plotfile"          : cfg.obs.basename + "_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "phaseuvdist.png",
            "overwrite"         :   True,
        },
        input=INPUT, output=OUTPUT,
        label="plot_phaseuvdist:: Diagnostic plot of phase with uvdist")

    # Diagnostic: amplitude vs phase of bp calibrator per antenna
    recipe.add("cab/casa_plotms", "plot_amp_v_phase_of_bp_calibrator",
        {
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
            "plotfile"          : cfg.obs.basename + "_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "phaseball.png",
            "overwrite"         :   True,
        },
        input=INPUT, output=OUTPUT,
        label="plot_phaseball:: Diagnostic plot of phaseball")

    # Diagnostic: amplitude vs frequency of bp calibrator
    recipe.add("cab/casa_plotms", "plot_amp_v_freq_of_bp_calibrator",
        {
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
            "plotfile"          : cfg.obs.basename + "_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "band.png",
            "overwrite"         :   True,
        },
        input=INPUT, output=OUTPUT,
        label="plot_amp_freq:: Diagnostic plot of band")

    # Diagnostic: phase vs time of bp calibrator
    # @oms points out slopes in this will indicate problems with
    # digitizer reference timing / probably also any uncorrected
    # antenna positions
    recipe.add("cab/casa_plotms", "plot_phase_vs_time_of_bp_calibrator",
        {
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
            "plotfile"          : cfg.obs.basename + "_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "phasevtime.png",
            "overwrite"         :   True,
        },
        input=INPUT, output=OUTPUT,
        label="plot_phase_time:: Diagnostic plot of phase with time")

    # Diagnostic: phase vs time of bp calibrator
    # For similar purposes as phase vs freq
    recipe.add("cab/casa_plotms", "plot_phase_vs_freq_of_bp_calibrator",
        {
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
            "plotfile"          : cfg.obs.basename + "_" +
                                  plot_name[bandpass_cal.name] + "_" +
                                  "phasevfreq.png",
            "overwrite"         :   True,
        },
        input=INPUT, output=OUTPUT,
        label="plot_phase_freq:: Diagnostic plot of phase with freq")

    # # this wastes disk space... who cares if the calibrators are in there or not
    # # it's actually very useful to keep them and inspect their solutions -
    # # DO NOT ADD THIS STEP
    # recipe.add("cab/casa_split", "split_calibrated_target_data",
    #     {
    #         "msname"        :   msname,
    #         "output_msname" :   msname[:-3]+"_deep2.ms",
    #         "field"         :   target,
    #     },
    # input=INPUT, output=OUTPUT,
    # label="split_target:: Split calibrated target data")

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
            },
            input=INPUT, output=OUTPUT,
            label="image_%d::wsclean" % target_field)

  # [Sphe] I don't think images of calibrators are needed. These are point sources, so images give us nothing we can't get from
  # gain/phase plots.

    # Diagnostic only: image bandpass
    recipe.add("cab/wsclean", "wsclean_bandpass",
        {
            "msname"            : cfg.obs.msfile,
            "column"            : cfg.wsclean_bandpass.column,
            "weight"            : "briggs %.2f"%(cfg.wsclean_bandpass.robust),
            "npix"              : cfg.obs.im_npix / 2, # don't need the full FOV
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
            "npix"              : cfg.obs.im_npix / 2, # don't need the full FOV
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

    # Add Flagging and 1GC steps
    onegcsteps = [
                   "convert",
                   "rfimask",
                   "autoflag",
                   "flag_bandstart",
                   "flag_bandend",
                   "flag_autocorrs",
                   "recompute_uvw",
                   "flag_baseline_phases",
                   "setjy",
                   "delaycal",
                   "phase0",
                   "bandpass",
                   "gaincal",
                   "fluxscale",
                   "applycal",
                   "autoflag_corrected_vis",
                   "plot_ampuvdist",
                   "plot_phaseuvdist",
                   "plot_phaseball",
                   "plot_amp_freq",
                   "plot_phase_time",
                   "plot_phase_freq",
                 ]
    onegcsteps += ["image_%d" % t for t in target_fields]
    onegcsteps += ["image_bandpass", "image_gain"] # diagnostic only

    # RUN FOREST RUN!!!
    # @simon The pipeline exceptions are now being handled in the recipe.run() function
    recipe.run(onegcsteps)

    ### Let the 2GC begin
    recipe = stimela.Recipe("2GC Pipeline", ms_dir=MSDIR)

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
                "thresh_pix"    : cfg.source_find.thresh_pix,
                "thresh_isl"    : cfg.source_find.thresh_isl,
                "port2tigger"   :   True,
            },
            input=INPUT, output=OUTPUT,
            label="source_find_%d:: Extract sources from previous round of cal" % target_field)

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
            label="stitch_cube_%d:: Stitch MFS image slices into a cube" % target_field)

        # Add SPIs
        recipe.add("cab/specfit", "add_SPIs_LSM0",
            {
                "image"                 :   "%s:output"%cubename,
                "output-spi-image"      :   "%s-spi.fits"%imname_prefix,
                "output-spi-error-image":   "%s-spi.error.fits"%imname_prefix,
                "input-skymodel"        :   "%s:output"%model_name, # model to which SPIs must be added
                "output-skymodel"       :   model_name,
                "tolerance-range"       :   cfg.specfit.tol,
                "freq0"                 :   cfg.obs.freq_0,    # reference frequency for SPI calculation
                "sigma-level"           :   cfg.specfit.sigma,
            },
            input=INPUT, output=OUTPUT,
            label="SPI_%d::Add SPIs to LSM" % target_field)

        # Selfcal and subtract brightest sources
        recipe.add("cab/calibrator", "Initial_Gjones_subtract_LSM0",
            {
                "skymodel"  :   "%s:output"%model_name,
                "label"     :   cfg.selfcal.label,
                "msname"    :   cfg.obs.msfile,
                "threads"   :   cfg.selfcal.ncpu,
                "column"    :   cfg.selfcal.column,
                "output-data"    :   cfg.selfcal.output,
                "Gjones"    :   cfg.selfcal.gjones,
                "Gjones-solution-intervals" : map(int, [float(cfg.obs.gain_sol_int[:-1]), cfg.obs.nchans / float(1000)]),
                "DDjones-smoothing-intervals" :  cfg.selfcal.ddjones_smoothing,
                # TODO: MeerKAT beams need to go in this section
                "Ejones"    :   cfg.selfcal.ejones,
                "beam-files-pattern" : cfg.selfcal.beam_files_pattern,
                "beam-l-axis" : cfg.selfcal.beam_l_axis,
                "beam-m-axis" : cfg.selfcal.beam_m_axis,
                "Gjones-ampl-clipping"  :   cfg.selfcal.gjones_ampl_clipping,
                "Gjones-ampl-clipping-low"  :   cfg.selfcal.gjones_ampl_clipping_low,
                "Gjones-ampl-clipping-high"  :   cfg.selfcal.gjones_ampl_clipping_high,
            },
            input=INPUT, output=OUTPUT,
            label="SELFCAL0_%d:: Calibrate and subtract LSM0" % target_field)

        #make another mfs image
        imname_prefix = cfg.obs.basename + "_SC0_" + plot_name[target.name]
        recipe.add("cab/wsclean", "wsclean_SC0_%d" % target_field,
            {
                "msname"            : cfg.obs.msfile,
                "column"            : cfg.wsclean_selfcal.column,
                "weight"            : "briggs %2.f"%(cfg.wsclean_selfcal.robust),
                "npix"              : cfg.obs.im_npix,
                "cellsize"          : cfg.obs.angular_resolution*cfg.obs.sampling,
                "clean_iterations"  : cfg.wsclean_selfcal.clean_iterations,
                "mgain"             : cfg.wsclean_selfcal.mgain,
                "channelsout"       : cfg.obs.im_numchans,
                "joinchannels"      : cfg.wsclean_selfcal.joinchannels,
                "field"             : str(target_field),
                "name"              : imname_prefix,
            },
            input=INPUT, output=OUTPUT,
            label="image_SC0_%d::wsclean" % target_field)

        #create a mask for this round of selfcal
        imname_mfs = imname_prefix + "-MFS-image.fits"
        maskname = imname_prefix + "_MASK"
        recipe.add("cab/cleanmask", "make_clean_mask",
            {
                "image"     :   '%s:output'%imname_mfs,
                "output"    :   maskname,
                "sigma"     :   cfg.cleanmask.sigma,
                "iters"     :   cfg.cleanmask.iters,
                "boxes"    :   cfg.cleanmask.kernel,
            },
            input=INPUT, output=OUTPUT,
            label="MSK_SC0_%d::Make clean mask" % target_field)

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
    twogcsteps = ["prepms",
                "move_corrdata_to_data",
               ]

    for target_field in target_fields:
        twogcsteps += ["source_find_%d" % target_field,
                       "stitch_cube_%d" % target_field,
                       "SPI_%d" % target_field,
                       "SELFCAL0_%d" % target_field,
                       "image_SC0_%d" % target_field,
                       "MSK_SC0_%d" % target_field,
                      ]

    # diagnostic only:
    twogcsteps += ["image_stokesv_residue_%d" % t for t in target_fields]

    recipe.run(twogcsteps)

