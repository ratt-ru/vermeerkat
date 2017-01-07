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
from functools import reduce

from vermeerkat.config import configuration
from vermeerkat.observation import (query_recent_observations,
    download_observation)

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
    bw_per_image_slice = 42.8e6 
    im_numchans = int(np.ceil(o["ChannelWidth"] * nchans / bw_per_image_slice))

    bandpass_cal_candidates = []
    gain_cal_candidates = []
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
        if scan[1] in source_name:
            continue
        source_name += [scan[1]]
        if "bpcal" in scan[2]:
            bandpass_cal_candidates += [si]
        elif "gaincal" in scan[2]:
            gain_cal_candidates += [si]
        elif "target" in scan[2]:
            targets += [si]
        else:
            vermeerkat.log.warn("Not using observed source %s" % source_string[0])

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
    target_coordinates = [scans[ti][3] for ti, t in enumerate(targets)]
    gaincal_coordinates = [scans[gi][3] for gi, g in enumerate(gain_cal_candidates)]
    mean_target_position = np.mean(np.array(target_coordinates), axis=0)
    lmdistances = np.sum((np.array(gaincal_coordinates) - mean_target_position)**2,
                         axis=0)
    gain_cal = np.argmin(lmdistances)
   
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

    vermeerkat.log.info("Will use '%s' (%d) as a gain calibrator" % 
        (source_name[gain_cal], gain_cal))
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
    vermeerkat.log.info("Will use firstpass aoflagger strategy: %s" % cfg.aoflagger.firstpass_strategy_file)
    vermeerkat.log.info("Will use secondpass aoflagger strategy: %s" % cfg.aoflagger.secondpass_strategy_file)
    vermeerkat.log.info("Will use RFI mask: %s" % cfg.rfimask.rfi_mask_file)
    vermeerkat.log.info("Correlator integration interval recorded as: %.2f secs" % correlator_integration_time)
    vermeerkat.log.info("MFS maps will contain %.3f Mhz per slice" % (bw_per_image_slice / 1e6))
    vermeerkat.log.info("Maps will cover %.3f square degrees at angular resolution %.3f asec" % 
        (fov / 3600.0, angular_resolution))
    vermeerkat.log.warn("Assuming maximum baseline is %.2f meters" % telescope_max_baseline) 
    vermeerkat.log.info("Will use '%s' as reference antenna" % refant)
    vermeerkat.log.info("Observed band covers %.2f MHz +/- %.2f MHz averaged into %d channels" % 
        (freq_0 / 1e6, (nchans // 2) * chan_bandwidth / 1e6, nchans))

    # All good to go fire up the pipeline
    recipe = stimela.Recipe("Imaging Pipeline", ms_dir=MSDIR)
    
    # Convert
    recipe.add("cab/h5toms", "h5toms",
        {
            'hdf5files'  : [h5file],
            'output-ms'  : msfile,
            'full_pol'   : True,
        },
        input=INPUT, output=OUTPUT,
        label="convert::h5toms")

    # RFI and bad channel flagging
    recipe.add("cab/rfimasker", "mask_stuff",
        {
            "msname" : msfile,
            "mask"   : cfg.rfimask.rfi_mask_file,
        },
        input=INPUT, output=OUTPUT,
        label="mask::maskms")

    recipe.add("cab/autoflagger", "auto_flag_rfi",
        {
            "msname"    : msfile,
            "column"    : "DATA",
            "strategy"  : cfg.aoflagger.firstpass_strategy_file,
        },
        input=INPUT, output=OUTPUT,
        label="autoflag:: Auto Flagging ms")
    
    recipe.add("cab/casa_flagdata", "flag_bad_start_channels",
        {
            "msname"    :   msfile,
            "mode"      :   "manual",
            "field"     :   '',
            "spw"       :   '0:0~109',
            "autocorr"  :   True,
        },
    	input=OUTPUT, output=OUTPUT,
        label="flag_bandstart:: Flag start of band")

    recipe.add("cab/casa_flagdata", "flag_bad_end_channels",
        {
            "msname"    :   msfile,
            "mode"      :   "manual",
            "field"     :   '',
            "spw"       :   '0:3860~4095',
            "autocorr"  :   True,
        },
    	input=OUTPUT, output=OUTPUT,
        label="flag_bandend:: Flag end of band")
    
    # 1GC Calibration
    recipe.add("cab/casa_setjy", "init_flux_scaling",
        {
            "msname"        :   msfile,
            "field"         :   str(bandpass_cal),
            # Relies on a custom casa docker image with
            # southern calibrators
            "standard"      :   'Perley-Butler 2013',
            #"fluxdensity"   :   24.5372,
            #"spix"          :   -1.02239,
            #"reffreq"       :   '900MHz',
            "usescratch"    :   False,
            "scalebychan"   :   True,
            "spw"           :   '',
        },
    	input=OUTPUT, output=OUTPUT,
        label="setjy:: Initial flux density scaling")

    recipe.add("cab/casa_gaincal", "init_phase_cal",
        {
            "msname"        :   msfile,
            "caltable"      :   phasecal_table,
            "field"         :   str(bandpass_cal),
            "refant"        :   refant,
            "calmode"       :   'p',
            "solint"        :   gain_sol_int,
            "minsnr"        :   3,
        },
    	input=OUTPUT, output=OUTPUT,
        label="phase0:: Initial phase calibration")


    recipe.add("cab/casa_bandpass", "bandpass_cal",
        {   
            "msname"        :   msfile, 
            "caltable"      :   bpasscal_table,
            "field"         :   str(bandpass_cal),
            "spw"           :   '',
            "refant"        :   refant,
            "combine"       :   'scan',
            "solint"        :   str(bpcal_sol_int) + "s",
            "bandtype"      :   'B',
            "minblperant"   :   1,
            "minsnr"        :   3,
            "gaintable"     :   [phasecal_table],
        },
    	input=OUTPUT, output=OUTPUT,
        label="bandpass:: First bandpass calibration")


    recipe.add("cab/casa_gaincal", "main_gain_calibration",
        {
            "msname"       :   msfile,
            "caltable"     :   ampcal_table,
            "field"        :   ",".join([str(x) for x in [bandpass_cal, gain_cal]]),
            "spw"          :   '',
            "solint"       :   gain_sol_int,
            "refant"       :   refant,
            "gaintype"     :   'G',
            "calmode"      :   'ap',
            "solnorm"      :   False,
            "gaintable"    :   [phasecal_table,
                                bpasscal_table],
            "interp"       :   ['linear','linear','nearest'],
        },
    	input=OUTPUT, output=OUTPUT,
        label="gaincal:: Gain calibration")


    recipe.add("cab/casa_fluxscale", "casa_fluxscale",
        {
            "msname"        :   msfile,
            "caltable"      :   ampcal_table,
            "fluxtable"     :   fluxcal_table,
            "reference"     :   [str(bandpass_cal)],
            "transfer"      :   [str(gain_cal)],
            "incremental"   :   False,
        },
    	input=OUTPUT, output=OUTPUT,
        label="fluxscale:: Setting Fluxscale")

    recipe.add("cab/casa_applycal", "apply_calibration", 
        {
            "msname"        :   msfile,
            "field"         :   ",".join([str(x) for x in (targets + [bandpass_cal, gain_cal])]),
            "gaintable"     :   [phasecal_table, bpasscal_table, fluxcal_table],
            "gainfield"     :   [bandpass_cal, bandpass_cal, gain_cal],
            "spwmap"        :   [[], [], []],
            "parang"        :   True,
        },
    	input=OUTPUT, output=OUTPUT,
        label="applycal:: Apply calibration solutions to target")

    recipe.add("cab/autoflagger", "auto_flag_rfi_corrected_vis",
        {
            "msname"    : msfile,
            "column"    : "CORRECTED_DATA",
            "strategy"  : cfg.aoflagger.secondpass_strategy_file,
        },
        input=INPUT, output=OUTPUT,
        label="autoflag_corrected_vis:: Auto Flagging calibrated visibilities")

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
        imname = basename + "_" + source_name[ti]
        recipe.add("cab/wsclean", "wsclean_%d" % ti,
            {
                "msname"            : msfile,
                "column"            : 'CORRECTED_DATA',
                "weight"            : 'briggs',
                "robust"            : 0,
                "npix"              : im_npix,
                "cellsize"          : angular_resolution,
                "clean_iterations"  : 1000,
                "mgain"             : 0.9,
                "channelsout"       : im_numchans,
                "joinchannels"      : True,
                "field"             : str(ti),
                "name"              : imname,
            },
            input=OUTPUT, output=OUTPUT,
            label="image_%d::wsclean" % ti)

    # Diagnostic only: image bandpass
    recipe.add("cab/wsclean", "wsclean_bandpass",
        {
            "msname"            : msfile,
            "column"            : 'CORRECTED_DATA',
            "weight"            : 'briggs',
            "robust"            : 0,
            "npix"              : im_npix / 5, # don't need the full FOV
            "cellsize"          : angular_resolution,
            "clean_iterations"  : 10000,
            "mgain"             : 0.9,
            "channelsout"       : im_numchans,
            "joinchannels"      : True,
            "field"             : str(bandpass_cal),
            "name"              : basename + "_bp_" + source_name[bandpass_cal],
        },
        input=OUTPUT, output=OUTPUT,
        label="image_bandpass::wsclean")

    # Diagnostic only: image gaincal
    recipe.add("cab/wsclean", "wsclean_gain",
        {
            "msname"            : msfile,
            "column"            : 'CORRECTED_DATA',
            "weight"            : 'briggs',
            "robust"            : 0,
            "npix"              : im_npix / 5, # don't need the full FOV
            "cellsize"          : angular_resolution,
            "clean_iterations"  : 10000,
            "mgain"             : 0.9,
            "channelsout"       : im_numchans,
            "joinchannels"      : True,
            "field"             : str(gain_cal),
            "name"              : basename + "_gc_" + source_name[gain_cal],
        },
        input=OUTPUT, output=OUTPUT,
        label="image_gain::wsclean")

    try:
        recipe.run("convert "
                   "mask "
                   "autoflag "
                   "flag_bandstart "
                   "flag_bandend "
                   "setjy "
                   "phase0 "
                   "bandpass "
                   "gaincal "
                   "fluxscale "
                   "applycal "
                   "autoflag_corrected_vis".split() + 
                   ["image_%d" % ti for ti in targets] +
                   ["image_bandpass", "image_gain"])

    except stimela.PipelineException as e:
        print 'completed {}'.format([c.label for c in e.completed])
        print 'failed {}'.format(e.failed.label)
        print 'remaining {}'.format([c.label for c in e.remaining])
        raise
