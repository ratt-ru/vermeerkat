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
    recipe = stimela.Recipe("Initial flagging Engine", ms_dir=MSDIR)

    # RFI and bad channel flagging
    recipe.add("cab/rfimasker", "mask_knownrfi",
               {
                   "msname": cfg.obs.msfile,
                   "mask": cfg.rfimask.rfi_mask_file,
               },
               input=INPUT, output=OUTPUT,
               label="rfimask::maskms")

    recipe.add("cab/autoflagger", "auto_flag_rfi",
               {
                   "msname": cfg.obs.msfile,
                   "column": cfg.autoflag.column,
                   "strategy": cfg.autoflag.strategy_file,
               },
               input=INPUT, output=OUTPUT,
               label="autoflag:: Auto Flagging ms")

    recipe.add("cab/casa_flagdata", "flag_bad_start_channels",
               {
                   "msname": cfg.obs.msfile,
                   "mode": cfg.flag_bandstart.mode,
                   "field": cfg.flag_bandstart.field,
                   "spw": cfg.flag_bandstart.spw,
                   "autocorr": cfg.flag_bandstart.autocorr,
               },
               input=INPUT, output=OUTPUT,
               label="flag_bandstart:: Flag start of band")

    recipe.add("cab/casa_flagdata", "flag_bad_end_channels",
               {
                   "msname": cfg.obs.msfile,
                   "mode": cfg.flag_bandend.mode,
                   "field": cfg.flag_bandend.field,
                   "spw": cfg.flag_bandend.spw,
                   "autocorr": cfg.flag_bandend.autocorr,
               },
               input=INPUT, output=OUTPUT,
               label="flag_bandend:: Flag end of band")

    recipe.add("cab/casa_flagdata", "flag_autocorrs",
               {
                   "msname": cfg.obs.msfile,
                   "mode": cfg.flag_autocorrs.mode,
                   "field": cfg.flag_autocorrs.field,
                   "spw": cfg.flag_autocorrs.spw,
                   "autocorr": cfg.flag_autocorrs.autocorr,
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