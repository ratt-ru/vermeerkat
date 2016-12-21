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

import stimela

INPUT = "input"
OUTPUT = "output"
MSDIR  = "msdir"

msname = "1466309464_1284_1024ch.full_pol.ms"
PREFIX = msname[:-3]

# Calibration tables
PHASECAL_TABLE = "${OUTPUT}/%s.G0"%PREFIX
AMPCAL_TABLE = "${OUTPUT}/%s.G1"%PREFIX
FLUXCAL_TABLE = "${OUTPUT}/%s.fluxscale"%PREFIX
BPASSCAL_TABLE = "${OUTPUT}/%s.B0"%PREFIX

rfi_mask_file = "rfi_mask.pickle"
strategy_file = "first_pass_meerkat_ar1.rfis"
bandpass_cal = "0408-65"
amp_cal = "0252-712" 
phase_cal = "2"
flux_cal = "0"
target = "3"
cal_model = "3C138_L.im"
gain_cal = "1"
refant = "m008"


# So that we can access GLOBALS pass through to the run command
stimela.register_globals()

recipe = stimela.Recipe("Imaging Pipeline", ms_dir=MSDIR)

recipe.add("cab/rfimasker", "mask_stuff",
    {
        "msname" : msname,
        "mask"   : rfi_mask_file,
    },
    input=INPUT, output=OUTPUT,
    label="mask::maskms")

recipe.add("cab/autoflagger", "auto_flag_rfi",
    {
        "msname"    : msname,
        "column"    : "DATA",
        "field"     : "0,1,2",
        "strategy"  : strategy_file,
    },
    input=INPUT, output=OUTPUT,
    label="autoflag:: Auto Flagging ms")

recipe.add("cab/casa_flagdata", "flag_bad_channels",
    {
        "msname"    :   msname,
        "mode"      :   "manual",
        "field"     :   '',
        "spw"       :   '0:0~109',
        "autocorr"  :   True,
    },
    input=INPUT, output=OUTPUT,
    label="flag_band:: Flag start of band")

## 1GC Calibration
recipe.add("cab/casa_setjy", "init_flux_scaling",
    {
        "msname"        :   msname,
        "field"         :   bandpass_cal,
        "standard"      :   'manual',
        "fluxdensity"   :   24.5372,
        "spix"          :   -1.02239,
        "reffreq"        :   '900MHz',
        "usescratch"    :   False,
        "scalebychan"   :   True,
        "spw"           :   '',
    },
    input=INPUT, output=OUTPUT,
    label="setjy:: Initial flux density scaling")

recipe.add("cab/casa_gaincal", "init_phase_cal",
    {
        "msname"        :   msname,
        "caltable"      :   PHASECAL_TABLE,
        "field"         :   bandpass_cal,
        "refant"        :   refant,
        "calmode"       :   'p',
        "solint"        :   '12s',
        "minsnr"        :   3,
    },
    input=INPUT, output=OUTPUT,
    label="phase0:: Initial phase calibration")


recipe.add("cab/casa_bandpass", "bandpass_cal",
    {   
        "msname"        :   msname, 
        "caltable"      :   BPASSCAL_TABLE,
        "field"         :   bandpass_cal,
        "spw"           :   '',
        "refant"        :   refant,
        "combine"       :   'scan',
        "solint"        :   '5min',
        "bandtype"      :   'B',
        "minblperant"   :   1,
        "gaintable"     :   [PHASECAL_TABLE],
    },
    input=INPUT, output=OUTPUT,
    label="bandpass:: First bandpass calibration")


recipe.add("cab/casa_gaincal", "main_gain_calibration",
    {
        "msname"        :   msname,
        "caltable"      :   AMPCAL_TABLE,
         "field"        :   "%s,%s"%(bandpass_cal, amp_cal),
         "spw"          :   '',
         "solint"       :   'inf',
         "refant"       :   refant,
         "gaintype"     :   'G',
         "calmode"      :   'ap',
         "solnorm"      :   False,
         "gaintable"    :   [PHASECAL_TABLE,
                             BPASSCAL_TABLE],
         "interp"       :   ['linear','linear','nearest'],
    },
    input=INPUT, output=OUTPUT,
    label="gaincal:: Gain calibration")


recipe.add("cab/casa_fluxscale", "casa_fluxscale",
    {
        "msname"        :   msname,
        "caltable"      :   AMPCAL_TABLE,
        "fluxtable"     :   FLUXCAL_TABLE,
        "reference"     :   [bandpass_cal],
        "transfer"      :   [amp_cal],
        "incremental"   :   False,
    },
        input=INPUT, output=OUTPUT,
        label="fluxscale:: Setting Fluxscale")

recipe.add("cab/casa_applycal", "apply_calibration", 
    {
        "msname"        :   msname,
        "field"         :   target,
        "gaintable"     :   [PHASECAL_TABLE, BPASSCAL_TABLE, FLUXCAL_TABLE],
        "gainfield"     :   [bandpass_cal, bandpass_cal, amp_cal],
        "spwmap"        :   [[], [], []],
        "parang"        :   True,
    },
    input=INPUT, output=OUTPUT,
    label="applycal:: Apply calibration solutions to target")

recipe.add("cab/casa_split", "split_calibrated_target_data",
    {
        "msname"        :   msname,
        "output_msname" :   msname[:-3]+"_deep2.ms",
        "field"         :   target,
    },
    input=INPUT, output=OUTPUT,
    label="split_target:: Split calibrated target data")


try:
    recipe.run("flag_band setjy phase0 bandpass gaincal fluxscale applycal".split())
except stimela.PipelineException as e:
    print 'completed {}'.format([c.label for c in e.completed])
    print 'failed {}'.format(e.failed.label)
    print 'remaining {}'.format([c.label for c in e.remaining])
    raise
