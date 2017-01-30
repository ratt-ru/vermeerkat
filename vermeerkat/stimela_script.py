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
import numpy as np
import re
import stimela
import vermeerkat
import katdal
import glob
from vermeerkat.caltable_parser import read_caltable
from vermeerkat.caltable_parser import convert_pb_to_casaspi
from functools import reduce

from vermeerkat.config import configuration
from vermeerkat.observation import (query_recent_observations,
    download_observation)

# There isn't a Southern standard in CASA
# so construct a little database of them for reference
vermeerkat.log.info("Parsing calibrator table")
ref_table = os.path.dirname(
    os.path.abspath(vermeerkat.__file__)) + "/southern_calibrators.txt"
calibrator_db = read_caltable(ref_table)

vermeerkat.log.info("Found the following reference calibrators (in GHz format):")
for key in calibrator_db:
    name = key
    epoch = calibrator_db[name]["epoch"]
    ra = calibrator_db[name]["ra"]
    decl = calibrator_db[name]["decl"]
    ag = calibrator_db[name]["a_ghz"]
    bg = calibrator_db[name]["b_ghz"]
    cg = calibrator_db[name]["c_ghz"]
    dg = calibrator_db[name]["d_ghz"]
    vermeerkat.log.info("\t%s\tEpoch:%d\tRA:%3.2f\tDEC:%3.2f\t"
                       "a:%.4f\tb:%.4f\tc:%.4f\td:%.4f" %
        (name, epoch, ra, decl, ag, bg, cg, dg))

# Register directories

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

    # Calibrator tables
    delaycal_table = "%s.K0" % (basename)
    phasecal_table = "%s.G0" % (basename)
    ampcal_table = "%s.G1" % (basename)
    fluxcal_table = "%s.fluxscale" % (basename)
    bpasscal_table = "%s.B0" % (basename)

    #Observation properties
    refant = o["RefAntenna"]
    correlator_integration_time = o["DumpPeriod"]
    gain_sol_int = str(correlator_integration_time * 3) + "s"
    freq_0 = o["CenterFrequency"]
    nchans = o["NumFreqChannels"]
    chan_bandwidth = o["ChannelWidth"]
    lambda_max = (299792458.0 / (freq_0 + (nchans // 2) * chan_bandwidth))
    telescope_max_baseline = 8e4
    angular_resolution = np.rad2deg(lambda_max /
                                    telescope_max_baseline *
                                    1.220) * 3600 #nyquest rate in arcsecs
    fov = 3 * 3600 # 3 degrees (ample to cover MeerKAT primary beam in L-Band)
    im_npix = fov / angular_resolution
    bw_per_image_slice = 100.0e6
    im_numchans = int(np.ceil(o["ChannelWidth"] * nchans / bw_per_image_slice))

    bandpass_cal_candidates = []
    gain_cal_candidates = []
    bandpass_coordinates = []
    gaincal_coordinates = []
    target_coordinates = []
    targets = []
    source_name = []

    d = katdal.open(INPUT + "/" + h5file)
    scans = [(scan,
              target.name,
              target.tags,
              target.radec(),
              d.timestamps[-1] - d.timestamps[0])
             for (scan,
                  state,
                  target) in d.scans() if state == "track"]

    for si, scan in enumerate(scans):
        catagorized_source = False
        if scan[1] in source_name:
            continue
        if "bpcal" in scan[2]:
            bandpass_cal_candidates += [len(source_name)]
            catagorized_source = True
        if "gaincal" in scan[2]:
            gain_cal_candidates += [len(source_name)]
            catagorized_source = True
        if "target" in scan[2]:
            targets += [len(source_name)]
            catagorized_source = True
        if not catagorized_source:
            vermeerkat.log.warn("Not using observed source %s" % source_string[0])
            continue
        source_name += [scan[1]]

    if len(targets) < 1:
        raise RuntimeError("Observation %s does not have any "
                           "targets" % o["ProductName"])
    if len(gain_cal_candidates) < 1:
        raise RuntimeError("Observation %s does not have any "
                           "gain calibrators" % o["ProductName"])
    if len(bandpass_cal_candidates) < 1:
        raise RuntimeError("Observation %s does not have any "
                           "bandpass calibrators" % o["ProductName"])

    # select gaincal candidate closest to the centre of the cluster of targets
    target_coordinates = [[s[3] for s in scans if s[1] == source_name[t]][0] for t in targets]
    gaincal_coordinates = [[s[3] for s in scans if s[1] == source_name[t]][0] for t in gain_cal_candidates]
    mean_target_position = np.mean(np.array(target_coordinates), 
                                   axis=0)

    lmdistances = np.sum((np.array(gaincal_coordinates) -
                         mean_target_position)**2,
                         axis=1)

    gain_cal = gain_cal_candidates[np.argmin(lmdistances)]
    # never use the gain calibrator for a bandpass calibrator
    bandpass_cal_candidates = [x for x in bandpass_cal_candidates if x != gain_cal]
    bp_cand_obs_len = [reduce((lambda x, y: x + y),
                              list(map((lambda q: q[4]),
                                   filter((lambda z: z[1] == source_name[cand]),
                                          scans))))
                       for cand in bandpass_cal_candidates]
    bandpass_cal = bandpass_cal_candidates[np.argmax(bp_cand_obs_len)]
    bandpass_scan_times = [min(list(map((lambda q: q[4]),
                               filter((lambda z: z[1] == source_name[cand]),
                                      scans))))
                           for cand in bandpass_cal_candidates]
    bpcal_sol_int = bandpass_scan_times[bandpass_cal_candidates.index(bandpass_cal)]

    # compute the observation time spent on gain calibrator:
    gaincal_scan_times = [reduce((lambda x, y: x + y),
                                list(map((lambda q: q[4]),
                                     filter((lambda z: z[1] == source_name[gc]),
                                            scans))))
                         for gc in gain_cal_candidates]
 

    # compute the observation time spent on targets:
    target_scan_times = [reduce((lambda x, y: x + y),
                                list(map((lambda q: q[4]),
                                     filter((lambda z: z[1] == source_name[t]),
                                            scans))))
                         for t in targets]


    # print out some useful statistics:
    vermeerkat.log.info("The following targets were observed:")
    for ti, t in enumerate(targets):
        obs_hr = int(np.floor(target_scan_times[ti] / 3600.0))
        obs_min = int(np.floor((target_scan_times[ti] - obs_hr * 3600.0) / 60.0))
        obs_sec = target_scan_times[ti] - obs_hr * 3600.0 - obs_min * 60.0
        vermeerkat.log.info("\t %s - total observation time %d h %d mins %.2f secs" %
            (source_name[t],
             obs_hr,
             obs_min,
             obs_sec))

    vermeerkat.log.info("Will use '%s' (%d) as a closest gain calibrator (total "
        "observation time: %.2f mins)" %
        (source_name[gain_cal], 
         gain_cal, 
         gaincal_scan_times[gain_cal_candidates.index(gain_cal)] / 60.0))
    vermeerkat.log.info("Will use '%s' (%d) as a bandpass calibrator (total "
        "observation time: %.2f mins, minimum scan time: %.2f mins)" %
        (source_name[bandpass_cal],
         bandpass_cal,
         bp_cand_obs_len[bandpass_cal_candidates.index(bandpass_cal)] / 60.0,
         bpcal_sol_int / 60.0))
    vermeerkat.log.info("Will image the following targets: '%s' (%s)" %
        (",".join([source_name[t] for t in targets]),
         ",".join([str(t) for t in targets])))
    vermeerkat.log.info("Will write ms file to %s" % msfile)
    vermeerkat.log.info("Will use firstpass aoflagger strategy: %s" % cfg.autoflag.strategy_file)
    vermeerkat.log.info("Will use secondpass aoflagger strategy: %s" % cfg.autoflag_corrected_vis.strategy_file)
    vermeerkat.log.info("Will use RFI mask: %s" % cfg.rfimask.rfi_mask_file)
    vermeerkat.log.info("Correlator integration interval recorded as: %.2f secs" % correlator_integration_time)
    vermeerkat.log.info("MFS maps will contain %.3f MHz per slice" % (bw_per_image_slice / 1e6))
    vermeerkat.log.info("Maps will cover %.3f square degrees at angular resolution %.3f asec" %
        (fov / 3600.0, angular_resolution))
    vermeerkat.log.warn("Assuming maximum baseline is %.2f meters" % telescope_max_baseline)
    vermeerkat.log.info("Will use '%s' as reference antenna" % refant)
    vermeerkat.log.info("Observed band covers %.2f MHz +/- %.2f MHz averaged into %d channels" %
        (freq_0 / 1e6, (nchans // 2) * chan_bandwidth / 1e6, nchans))

    # All good to go fire up the pipeline
    recipe = stimela.Recipe("1GC Pipeline", ms_dir=MSDIR)

    # Convert
    recipe.add("cab/h5toms", "h5toms",
        {
            'hdf5files'  : [h5file],
            'output-ms'  : msfile,
            'full_pol'   : cfg.h5toms.full_pol,
        },
        input=INPUT, output=OUTPUT,
        label="convert::h5toms")

    recipe.add("cab/msutils", "msutils",
        {
            'command'    : 'prep',
            'msname'     : msfile,
        },
        input=OUTPUT, output=OUTPUT,
        label="prepms::Adds flagsets")

    # RFI and bad channel flagging
    recipe.add("cab/rfimasker", "mask_stuff",
        {
            "msname" : msfile,
            "mask"   : cfg.rfimask.rfi_mask_file,
        },
        input=INPUT, output=OUTPUT,
        label="rfimask::maskms")

    recipe.add("cab/autoflagger", "auto_flag_rfi",
        {
            "msname"    : msfile,
            "column"    : cfg.autoflag.column,
            "strategy"  : cfg.autoflag.strategy_file,
        },
        input=INPUT, output=OUTPUT,
        label="autoflag:: Auto Flagging ms")

    # Custom user specified flags
    # Casa does not sensibly do nothing if the defaults are specified
    # removing this step
#    userflags = {"msname" : msfile}
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
#        input=OUTPUT, output=OUTPUT,
#        label="flag_userflags:: Specify any additional flags from user")

    recipe.add("cab/casa_flagdata", "flag_bad_start_channels",
        {
            "msname"    :   msfile,
            "mode"      :   cfg.flag_bandstart.mode,
            "field"     :   cfg.flag_bandstart.field,
            "spw"       :   cfg.flag_bandstart.spw,
            "autocorr"  :   cfg.flag_bandstart.autocorr,
        },
        input=OUTPUT, output=OUTPUT,
        label="flag_bandstart:: Flag start of band")

    recipe.add("cab/casa_flagdata", "flag_bad_end_channels",
        {
            "msname"    :   msfile,
            "mode"      :   cfg.flag_bandend.mode,
            "field"     :   cfg.flag_bandend.field,
            "spw"       :   cfg.flag_bandend.spw,
            "autocorr"  :   cfg.flag_bandend.autocorr,
        },
        input=OUTPUT, output=OUTPUT,
        label="flag_bandend:: Flag end of band")

    recipe.add("cab/casa_flagdata", "flag_autocorrs",
        {
            "msname"    :   msfile,
            "mode"      :   cfg.flag_autocorrs.mode,
            "field"     :   cfg.flag_autocorrs.field,
            "spw"       :   cfg.flag_autocorrs.spw,
            "autocorr"  :   cfg.flag_autocorrs.autocorr,
        },
        input=OUTPUT, output=OUTPUT,
        label="flag_autocorrs:: Flag auto correlations")

    # @mauch points out some MeerKAT h5 files contain flipped uvw
    # coordinates. Better to recalculate this from ANTENNAS
    # and be dead sure things will always work!
    recipe.add("cab/casa_fixvis", "fixvis",
        {
            "vis"       :   msfile,
            "outputvis" :   msfile, # into same ms please
            "reuse"     :   cfg.casa_fixvis.reuse,
        },
        input=OUTPUT, output=OUTPUT,
        label="recompute_uvw:: Recompute MeerKAT uvw coordinates")

    #Run SetJY with our database of southern calibrators
    if source_name[bandpass_cal] not in calibrator_db:
        raise RuntimeError("Looks like your flux reference '%s' is not "
                           "in our standard. Try pulling the latest "
                           "VermeerKAT or if you have done "
                           "so report this issue" % source_name[bandpass_cal])

    aghz = calibrator_db[source_name[bandpass_cal]]["a_ghz"]
    bghz = calibrator_db[source_name[bandpass_cal]]["b_ghz"]
    cghz = calibrator_db[source_name[bandpass_cal]]["c_ghz"]
    dghz = calibrator_db[source_name[bandpass_cal]]["d_ghz"]

    # Find the brightness at reference frequency
    I, a, b, c, d = convert_pb_to_casaspi(freq_0 / 1e9 -
                                          (nchans // 2) * chan_bandwidth / 1e9,
                                          freq_0 / 1e9 +
                                          (nchans // 2) * chan_bandwidth / 1e9,
                                          freq_0 / 1e9,
                                          aghz,
                                          bghz,
                                          cghz,
                                          dghz)

    vermeerkat.log.info("Using bandpass calibrator %s "
                        "with brightness of %.4f Jy "
                        "(spix = [%.6f, %.6f, %.6f, %.6f]) "
                        "(@ %.2f MHz) "
                        "as the flux scale reference" %
                        (source_name[bandpass_cal],
                         I, a, b, c, d,
                         freq_0 / 1e6))
    # 1GC Calibration
    recipe.add("cab/casa_setjy", "init_flux_scaling",
        {
            "msname"        :   msfile,
            "field"         :   str(bandpass_cal),
            "standard"      :   cfg.setjy.standard,
            "fluxdensity"   :   I,
            "spix"          :   [a, b, c, d],
            "reffreq"       :   "%.2fGHz" % (freq_0 / 1e9),
            "usescratch"    :   cfg.setjy.usescratch,
            "scalebychan"   :   cfg.setjy.scalebychan,
            "spw"           :   cfg.setjy.spw,
        },
        input=OUTPUT, output=OUTPUT,
        label="setjy:: Initial flux density scaling")

    # @mauch points out some antenna positions may
    # be slightly off. This will cause a noticible
    # decorrelation on the longest baselines, so
    # better calibrate for this one (we can use the
    # bright bandpass and gain calibrator sources).
    recipe.add("cab/casa_gaincal", "delay_cal",
        {
            "msname"        : msfile,
            "field"         : str(gain_cal),
            "gaintype"      : cfg.delaycal.gaintype,
            "solint"        : cfg.delaycal.solint,
            "minsnr"        : cfg.delaycal.minsnr,
            "refant"        : refant,
            "caltable"      : delaycal_table,
            "calmode"       : cfg.delaycal.calmode,
        },
        input=OUTPUT, output=OUTPUT,
        label="delaycal:: Delay calibration")

    # The bandpass calibrator may vary in phase over time
    # so lets do a preliminary phase cal to correct for this
    recipe.add("cab/casa_gaincal", "init_phase_cal",
        {
            "msname"        :   msfile,
            "caltable"      :   phasecal_table,
            "field"         :   str(bandpass_cal),
            "refant"        :   refant,
            "calmode"       :   cfg.phase0.calmode,
            "solint"        :   gain_sol_int,
            "solnorm"       :   cfg.phase0.solnorm,
            "minsnr"        :   cfg.phase0.minsnr,
            "gaintable"     :   [delaycal_table],
        },
        input=OUTPUT, output=OUTPUT,
        label="phase0:: Initial phase calibration")

    # Then a bandpass calibration, lets work out a solution
    # per scan interval to see if things vary signifcantly
    # over time... this will be very bad for your reduction.
    recipe.add("cab/casa_bandpass", "bandpass_cal",
        {
            "msname"        :   msfile,
            "caltable"      :   bpasscal_table,
            "field"         :   str(bandpass_cal),
            "spw"           :   cfg.bandpass.spw,
            "refant"        :   refant,
            "solnorm"       :   cfg.bandpass.solnorm,
            "combine"       :   cfg.bandpass.combine,
            "solint"        :   str(bpcal_sol_int) + "s",
            "bandtype"      :   cfg.bandpass.bandtype,
            "minblperant"   :   cfg.bandpass.minblperant,
            "minsnr"        :   cfg.bandpass.minsnr,
            "gaintable"     :   [delaycal_table,
                                 phasecal_table],
            "interp"        :   cfg.bandpass.interp,
        },
        input=OUTPUT, output=OUTPUT,
        label="bandpass:: Bandpass calibration")

    # Finally we do a second order correction on the gain
    # cal source that is closest to the target fields
    recipe.add("cab/casa_gaincal", "main_gain_calibration",
        {
            "msname"       :   msfile,
            "caltable"     :   ampcal_table,
            "field"        :   ",".join([str(x) for x in [bandpass_cal, gain_cal]]),
            "spw"          :   cfg.gaincal.spw,
            "solint"       :   gain_sol_int,
            "refant"       :   refant,
            "gaintype"     :   cfg.gaincal.gaintype,
            "calmode"      :   cfg.gaincal.calmode,
            "solnorm"      :   cfg.gaincal.solnorm,
            "gaintable"    :   [delaycal_table,
                                phasecal_table,
                                bpasscal_table],
            "interp"       :   cfg.gaincal.interp,
        },
        input=OUTPUT, output=OUTPUT,
        label="gaincal:: Gain calibration")

    # Scale the gaincal solutions amplitude to that of the
    # bandpass calibrator (the flux scale reference)
    recipe.add("cab/casa_fluxscale", "casa_fluxscale",
        {
            "msname"        :   msfile,
            "caltable"      :   ampcal_table,
            "fluxtable"     :   fluxcal_table,
            "reference"     :   [str(bandpass_cal)],
            "transfer"      :   [str(gain_cal)],
            "incremental"   :   cfg.fluxscale.incremental,
        },
        input=OUTPUT, output=OUTPUT,
        label="fluxscale:: Setting Fluxscale")

    # Apply gain solutions to all fields including
    # the calibrators so that we can diagnose problems more
    # easily. Later steps depend on this so don't remove
    # the gain application to calibrators
    recipe.add("cab/casa_applycal", "apply_calibration",
        {
            "msname"        :   msfile,
            "field"         :   ",".join([str(x) for x in (targets + [bandpass_cal, gain_cal])]),
            "gaintable"     :   [delaycal_table,
                                 phasecal_table,
                                 bpasscal_table,
                                 fluxcal_table],
            "gainfield"     :   [gain_cal,
                                 bandpass_cal,
                                 bandpass_cal,
                                 gain_cal],
            "interp"        :   cfg.applycal.interp,
            "spwmap"        :   cfg.applycal.spwmap,
            "parang"        :   cfg.applycal.parang,
        },
        input=OUTPUT, output=OUTPUT,
        label="applycal:: Apply calibration solutions to target")

    # Lets try squash the last of the RFI with the autoflagger
    # TODO: the strategy used here may still need some fine tuning
    # Think we should make it a bit more aggressive
    recipe.add("cab/autoflagger", "auto_flag_rfi_corrected_vis",
        {
            "msname"    : msfile,
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
            "msname"                 : msfile,
            "data_column"            : cfg.flag_baseline_phases.data_column,
            "cal_field"              : str(bandpass_cal),
            "valid_phase_range"      : cfg.flag_baseline_phases.valid_phase_range,
            "max_invalid_datapoints" : cfg.flag_baseline_phases.max_invalid_datapoints,
            "output_dir"             : os.environ["HOME"] + "/" + OUTPUT + "/", #where is this thing mapped to inside the container
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
            "msname"            : msfile,
            "xaxis"             : cfg.plot_ampuvdist.xaxis,
            "yaxis"             : cfg.plot_ampuvdist.yaxis,
            "xdatacolumn"       : cfg.plot_ampuvdist.xdatacolumn,
            "ydatacolumn"       : cfg.plot_ampuvdist.ydatacolumn,
            "field"             : str(bandpass_cal),
            "correlation"       : cfg.plot_ampuvdist.correlation,
            "iteraxis"          : cfg.plot_ampuvdist.iteraxis,
            "avgchannel"        : cfg.plot_ampuvdist.avgchannel,
            "avgtime"           : cfg.plot_ampuvdist.avgtime,
            "coloraxis"         : cfg.plot_ampuvdist.coloraxis,
            "expformat"         : cfg.plot_ampuvdist.expformat,
            "exprange"          : cfg.plot_ampuvdist.exprange,
            "plotfile"          : basename + "_" +
                                  source_name[bandpass_cal] + "_" +
                                  "ampuvdist.png",
        },
        input=OUTPUT, output=OUTPUT,
        label="plot_ampuvdist:: Diagnostic plot of amplitude with uvdist")

    #Diagnostic: phase vs uv dist of the bandpass calibrator
    recipe.add("cab/casa_plotms", "plot_phase_v_uv_dist_of_bp_calibrator",
        {
            "msname"            : msfile,
            "xaxis"             : cfg.plot_phaseuvdist.xaxis,
            "yaxis"             : cfg.plot_phaseuvdist.yaxis,
            "xdatacolumn"       : cfg.plot_phaseuvdist.xdatacolumn,
            "ydatacolumn"       : cfg.plot_phaseuvdist.ydatacolumn,
            "field"             : str(bandpass_cal),
            "correlation"       : cfg.plot_phaseuvdist.correlation,
            "iteraxis"          : cfg.plot_phaseuvdist.iteraxis,
            "avgchannel"        : cfg.plot_phaseuvdist.avgchannel,
            "avgtime"           : cfg.plot_phaseuvdist.avgtime,
            "coloraxis"         : cfg.plot_phaseuvdist.coloraxis,
            "expformat"         : cfg.plot_phaseuvdist.expformat,
            "exprange"          : cfg.plot_phaseuvdist.exprange,
            "plotfile"          : basename + "_" +
                                  source_name[bandpass_cal] + "_" +
                                  "phaseuvdist.png",
        },
        input=OUTPUT, output=OUTPUT,
        label="plot_phaseuvdist:: Diagnostic plot of phase with uvdist")

    # Diagnostic: amplitude vs phase of bp calibrator per antenna
    recipe.add("cab/casa_plotms", "plot_amp_v_phase_of_bp_calibrator",
        {
            "msname"            : msfile,
            "xaxis"             : cfg.plot_phaseball.xaxis,
            "yaxis"             : cfg.plot_phaseball.yaxis,
            "xdatacolumn"       : cfg.plot_phaseball.xdatacolumn,
            "ydatacolumn"       : cfg.plot_phaseball.ydatacolumn,
            "field"             : str(bandpass_cal),
            "correlation"       : cfg.plot_phaseball.correlation,
            "iteraxis"          : cfg.plot_phaseball.iteraxis,
            "avgchannel"        : cfg.plot_phaseball.avgchannel,
            "avgtime"           : cfg.plot_phaseball.avgtime,
            "coloraxis"         : cfg.plot_phaseball.coloraxis,
            "expformat"         : cfg.plot_phaseball.expformat,
            "exprange"          : cfg.plot_phaseball.exprange,
            "plotfile"          : basename + "_" +
                                  source_name[bandpass_cal] + "_" +
                                  "phaseball.png",
        },
        input=OUTPUT, output=OUTPUT,
        label="plot_phaseball:: Diagnostic plot of phaseball")

    # Diagnostic: amplitude vs frequency of bp calibrator
    recipe.add("cab/casa_plotms", "plot_amp_v_freq_of_bp_calibrator",
        {
            "msname"            : msfile,
            "xaxis"             : cfg.plot_amp_freq.xaxis,
            "yaxis"             : cfg.plot_amp_freq.yaxis,
            "xdatacolumn"       : cfg.plot_amp_freq.xdatacolumn,
            "ydatacolumn"       : cfg.plot_amp_freq.ydatacolumn,
            "field"             : str(bandpass_cal),
            "correlation"       : cfg.plot_amp_freq.correlation,
            "iteraxis"          : cfg.plot_amp_freq.iteraxis,
            "avgchannel"        : cfg.plot_amp_freq.avgchannel,
            "avgtime"           : cfg.plot_amp_freq.avgtime,
            "coloraxis"         : cfg.plot_amp_freq.coloraxis,
            "expformat"         : cfg.plot_amp_freq.expformat,
            "exprange"          : cfg.plot_amp_freq.exprange,
            "plotfile"          : basename + "_" +
                                  source_name[bandpass_cal] + "_" +
                                  "band.png",
        },
        input=OUTPUT, output=OUTPUT,
        label="plot_amp_freq:: Diagnostic plot of band")

    # Diagnostic: phase vs time of bp calibrator
    # @oms points out slopes in this will indicate problems with
    # digitizer reference timing / probably also any uncorrected
    # antenna positions
    recipe.add("cab/casa_plotms", "plot_phase_vs_time_of_bp_calibrator",
        {
            "msname"            : msfile,
            "xaxis"             : cfg.plot_phase_time.xaxis,
            "yaxis"             : cfg.plot_phase_time.yaxis,
            "xdatacolumn"       : cfg.plot_phase_time.xdatacolumn,
            "ydatacolumn"       : cfg.plot_phase_time.ydatacolumn,
            "field"             : str(bandpass_cal),
            "correlation"       : cfg.plot_phase_time.correlation,
            "iteraxis"          : cfg.plot_phase_time.iteraxis,
            "avgchannel"        : cfg.plot_phase_time.avgchannel,
            "avgtime"           : cfg.plot_phase_time.avgtime,
            "coloraxis"         : cfg.plot_phase_time.coloraxis,
            "expformat"         : cfg.plot_phase_time.expformat,
            "exprange"          : cfg.plot_phase_time.exprange,
            "plotfile"          : basename + "_" +
                                  source_name[bandpass_cal] + "_" +
                                  "phasevtime.png",
        },
        input=OUTPUT, output=OUTPUT,
        label="plot_phase_time:: Diagnostic plot of phase with time")

    # Diagnostic: phase vs time of bp calibrator
    # For similar purposes as phase vs freq
    recipe.add("cab/casa_plotms", "plot_phase_vs_freq_of_bp_calibrator",
        {
            "msname"            : msfile,
            "xaxis"             : cfg.plot_phase_freq.xaxis,
            "yaxis"             : cfg.plot_phase_freq.yaxis,
            "xdatacolumn"       : cfg.plot_phase_freq.xdatacolumn,
            "ydatacolumn"       : cfg.plot_phase_freq.ydatacolumn,
            "field"             : str(bandpass_cal),
            "correlation"       : cfg.plot_phase_freq.correlation,
            "iteraxis"          : cfg.plot_phase_freq.iteraxis,
            "avgchannel"        : cfg.plot_phase_freq.avgchannel,
            "avgtime"           : cfg.plot_phase_freq.avgtime,
            "coloraxis"         : cfg.plot_phase_freq.coloraxis,
            "expformat"         : cfg.plot_phase_freq.expformat,
            "exprange"          : cfg.plot_phase_freq.exprange,
            "plotfile"          : basename + "_" +
                                  source_name[bandpass_cal] + "_" +
                                  "phasevfreq.png",
        },
        input=OUTPUT, output=OUTPUT,
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
    for ti in targets:
        imname = basename + "_1GC_" + source_name[ti]
        recipe.add("cab/wsclean", "wsclean_%d" % ti,
            {
                "msname"            : msfile,
                "column"            : cfg.wsclean_image.column,
                "weight"            : cfg.wsclean_image.weight,
                "robust"            : cfg.wsclean_image.robust,
                "npix"              : im_npix,
                "cellsize"          : angular_resolution,
                "clean_iterations"  : cfg.wsclean_image.clean_iterations,
                "mgain"             : cfg.wsclean_image.mgain,
                "channelsout"       : im_numchans,
                "joinchannels"      : cfg.wsclean_image.joinchannels,
                "field"             : str(ti),
                "name"              : imname,
            },
            input=OUTPUT, output=OUTPUT,
            label="image_%d::wsclean" % ti)

    # Diagnostic only: image bandpass
    recipe.add("cab/wsclean", "wsclean_bandpass",
        {
            "msname"            : msfile,
            "column"            : cfg.wsclean_bandpass.column,
            "weight"            : cfg.wsclean_bandpass.weight,
            "robust"            : cfg.wsclean_bandpass.robust,
            "npix"              : im_npix / 5, # don't need the full FOV
            "cellsize"          : angular_resolution,
            "clean_iterations"  : cfg.wsclean_bandpass.clean_iterations,
            "mgain"             : cfg.wsclean_bandpass.mgain,
            "channelsout"       : im_numchans,
            "joinchannels"      : cfg.wsclean_bandpass.joinchannels,
            "field"             : str(bandpass_cal),
            "name"              : basename + "_bp_" + source_name[bandpass_cal],
        },
        input=OUTPUT, output=OUTPUT,
        label="image_bandpass::wsclean")

    # Diagnostic only: image gaincal
    recipe.add("cab/wsclean", "wsclean_gain",
        {
            "msname"            : msfile,
            "column"            : cfg.wsclean_gain.column,
            "weight"            : cfg.wsclean_gain.weight,
            "robust"            : cfg.wsclean_gain.robust,
            "npix"              : im_npix / 5, # don't need the full FOV
            "cellsize"          : angular_resolution,
            "clean_iterations"  : cfg.wsclean_gain.clean_iterations,
            "mgain"             : cfg.wsclean_gain.mgain,
            "channelsout"       : im_numchans,
            "joinchannels"      : cfg.wsclean_gain.joinchannels,
            "field"             : str(gain_cal),
            "name"              : basename + "_gc_" + source_name[gain_cal],
        },
        input=OUTPUT, output=OUTPUT,
        label="image_gain::wsclean")

    try:
        # Add Flagging and 1GC steps
        onegcsteps = ["convert",
                      "prepms",
                      "rfimask",
                      "autoflag",
                      "flag_bandstart",
                      "flag_bandend",
                      "flag_autocorrs",
                      "recompute_uvw",
                      "setjy",
                      "delaycal",
                      "phase0",
                      "bandpass",
                      "gaincal",
                      "fluxscale",
                      "applycal",
                      "autoflag_corrected_vis",
                      "flag_baseline_phases",
                      "plot_ampuvdist",
                      "plot_phaseuvdist",
                      "plot_phaseball",
                      "plot_amp_freq",
                      "plot_phase_time",
                      "plot_phase_freq",
                     ]
        onegcsteps += ["image_%d" % ti for ti in targets]
        onegcsteps += ["image_bandpass", "image_gain"] # diagnostic only
 
        # RUN FOREST RUN!!!
        recipe.run(onegcsteps)

    except stimela.PipelineException as e:
        print 'completed {}'.format([c.label for c in e.completed])
        print 'failed {}'.format(e.failed.label)
        print 'remaining {}'.format([c.label for c in e.remaining])
        raise
 
    recipe = stimela.Recipe("2GC Pipeline", ms_dir=MSDIR)
    # Copy CORRECTED_DATA to DATA, so we can start selfcal
    recipe.add("cab/msutils", "shift_columns",
        {
            "command"           : cfg.move_corrdata_to_data.command,
            "msname"            : msfile,
            "fromcol"           : cfg.move_corrdata_to_data.fromcol,
            "tocol"             : cfg.move_corrdata_to_data.tocol,
        },
        input=OUTPUT, output=OUTPUT,
        label="move_corrdata_to_data::msutils")

    # Mark current flags as legacy
    recipe.add("cab/flagms", "clear_flags",
        {
            "msname"        :  msfile,
            "flagged-any"   :  cfg.flagset_saveas_legacy.flagged_any,
            "flag"          :  cfg.flagset_saveas_legacy.flag,
        },
        input=INPUT,   output=OUTPUT,
        label="flagset_saveas_legacy:: Update legacy flags")

    # Initial selfcal loop
    for ti in targets:
        # Extract sources in mfs clean image to build initial sky model
        imname_prefix = basename + "_1GC_" + source_name[ti]
        imname_mfs = imname_prefix + "-MFS-image.fits"
        model_prefix = basename + "_LSM0_" + source_name[ti]
        model_name = model_prefix + ".lsm.html"
        recipe.add("cab/sourcery", "extract_sources_%d" % ti,
            {
                "imagename"     : imname_mfs,
                "prefix"        : model_prefix,
                "pybdsm"        : cfg.source_find.pybdsm,
                "thresh_pix"    : cfg.source_find.thresh_pix,
                "thresh_isl"    : cfg.source_find.thresh_isl,
            },
            input=OUTPUT, output=OUTPUT,
            label="source_find_%d:: Extract sources from previous round of cal" % ti)

        # Stitch wsclean channel images into a cube
        cubename = basename + "_1GC_" + source_name[ti] + "-CLEAN_cube.fits"
        recipe.add("cab/fitstool", "fitstool",
            {
                "image-names"       : [os.path.basename(img)
                                       for img in glob.glob("%s/%s-*-image.fits" %
                                                            (OUTPUT,
                                                             imname_prefix))],
                "output"            : cubename,
                "stack"             : "%s:3" % cubename,
            },
            input=OUTPUT, output=OUTPUT,
            label="stitch_cube_%d:: Stitch MFS image slices into a cube" % ti)

        # Add SPIs
        recipe.add("cab/specfit", "add_SPIs_LSM0",
            {
                "image"     :   cubename,
                "make_spi"  :   cfg.specfit.make_spi,
                "tol"       :   cfg.specfit.tol,
                "spi_image" :   cfg.specfit.spi_image,
                "add_spi"   :   cfg.specfit.add_spi,
                "skymodel"  :   model_name, # model to which SPIs must be added
                "freq0"     :   freq_0,    # reference frequency for SPI calculation
                "sigma"     :   cfg.specfit.sigma,
            },
            input=OUTPUT, output=OUTPUT,
            label="SPI_%d::Add SPIs to LSM" % ti)

        # Selfcal and subtract brightest sources
        recipe.add("cab/calibrator", "Initial_Gjones_subtract_LSM0",
            {
                "skymodel"  :   model_name,
                "reset"     :   cfg.selfcal.reset,
                "label"     :   cfg.selfcal.label,
                "msname"    :   msfile,
                "ncpu"      :   cfg.selfcal.ncpu,
                "column"    :   cfg.selfcal.column,
                "output"    :   cfg.selfcal.output,
                "Gjones"    :   cfg.selfcal.gjones,
                "Gjones_intervals" : [gain_sol_int, nchans / float(1000)],
                "DDjones_smoothing" :  cfg.selfcal.ddjones_smoothing,
                # TODO: MeerKAT beams need to go in this section
                "Ejones"    :   cfg.selfcal.ejones,
                "beam_files_pattern" : cfg.selfcal.beam_files_pattern,
                "beam_l_axis" : cfg.selfcal.beam_l_axis,
                "beam_m_axis" : cfg.selfcal.beam_m_axis,
                "gjones_ampl_clipping"  :   cfg.selfcal.gjones_ampl_clipping,
                "args"  :   cfg.selfcal.args,
            },
            input=OUTPUT, output=OUTPUT,
            label="SELFCAL0_%d:: Calibrate and subtract LSM0" % ti)

        #make another mfs image
        imname_prefix = basename + "_SC0_" + source_name[ti]
        recipe.add("cab/wsclean", "wsclean_SC0_%d" % ti,
            {
                "msname"            : msfile,
                "column"            : cfg.wsclean_selfcal.column,
                "weight"            : cfg.wsclean_selfcal.weight,
                "robust"            : cfg.wsclean_selfcal.robust,
                "npix"              : im_npix,
                "cellsize"          : angular_resolution,
                "clean_iterations"  : cfg.wsclean_selfcal.clean_iterations,
                "mgain"             : cfg.wsclean_selfcal.mgain,
                "channelsout"       : im_numchans,
                "joinchannels"      : cfg.wsclean_selfcal.joinchannels,
                "field"             : str(ti),
                "name"              : imname_prefix,
            },
            input=OUTPUT, output=OUTPUT,
            label="image_SC0_%d::wsclean" % ti)

        #create a mask for this round of selfcal
        imname_mfs = imname_prefix + "-MFS-image.fits"
        maskname = imname_prefix + "_MASK"
        recipe.add("cab/cleanmask", "make_clean_mask",
            {
                "image"     :   imname_mfs,
                "outname"   :   maskname,
                "sigma"     :   cfg.cleanmask.sigma,
                "iters"     :   cfg.cleanmask.iters,
                "kernel"    :   cfg.cleanmask.kernel,
            },
            input=INPUT, output=OUTPUT,
            label="MSK_SC0_%d::Make clean mask" % ti)

    for ti in targets:
        # Extract sources in mfs clean image to build initial sky model
        imname_prefix = basename + "_STOKES_V_RESIDUE_" + source_name[ti]

        # it  is common belief that the Unverse is free of
        # stokes V for the most part. So any structure left
        # in it is probably calibration artifacts
        recipe.add("cab/wsclean", "wsclean_stokes_v",
        {
            "msname"            : msfile,
            "column"            : cfg.wsclean_v_residue.column,
            "weight"            : cfg.wsclean_v_residue.weight,
            "robust"            : cfg.wsclean_v_residue.robust,
            "npix"              : im_npix,
            "cellsize"          : angular_resolution,
            "clean_iterations"  : cfg.wsclean_v_residue.clean_iterations,
            "mgain"             : cfg.wsclean_v_residue.mgain,
            "channelsout"       : im_numchans,
            "joinchannels"      : cfg.wsclean_v_residue.joinchannels,
            "field"             : str(ti),
            "name"              : imname_prefix,
            "pol"               : cfg.wsclean_v_residue.pol,
        },
        input=OUTPUT, output=OUTPUT,
        label="image_stokesv_residue_%d::wsclean image STOKES "
              "V as diagnostic" % ti)

    try:
        # Initial selfcal loop
        twogcsteps = ["move_corrdata_to_data",
                    "flagset_saveas_legacy",
                   ]

        for ti in targets:
            twogcsteps += ["source_find_%d" % ti,
                           "stitch_cube_%d" % ti,
                           "SPI_%d" % ti,
                           "SELFCAL0_%d" % ti,
                           "image_SC0_%d" % ti,
                           "MSK_SC0_%d" % ti,
                          ]

        # diagnostic only:
        twogcsteps += ["image_stokesv_residue_%d" % ti for ti in targets]

        # RUN FOREST RUN!!!
        recipe.run(twogcsteps)

    except stimela.PipelineException as e:
        print 'completed {}'.format([c.label for c in e.completed])
        print 'failed {}'.format(e.failed.label)
        print 'remaining {}'.format([c.label for c in e.remaining])
        raise

