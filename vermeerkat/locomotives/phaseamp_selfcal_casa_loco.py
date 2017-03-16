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
    plot_name = kwargs["plot_name"]
    targets = kwargs["targets"]
    target_fields = kwargs["target_fields"]
    recipe = stimela.Recipe("Phase amplitude 2GC Pipeline", ms_dir=MSDIR)
    # Initial selfcal loop
    for target_field, target in zip(target_fields, targets):
        gain_cal_opts = {
            "msname": cfg.obs.msfile,
            "caltable": cfg.obs.casa_SC1_gain_table,
            "field": str(target_field),
            "spw": cfg.casa_sc1_gaincal.spw,
            "solint": cfg.obs.gain_sol_int,
            "refant": cfg.obs.refant,
            "gaintype": cfg.casa_sc1_gaincal.gaintype,
            "calmode": cfg.casa_sc1_gaincal.calmode,
            "solnorm": cfg.casa_sc1_gaincal.solnorm,
        }
        apply_cal_opts = {
            "msname": cfg.obs.msfile,
            "field": str(target_field),
            "gaintable": [cfg.obs.casa_SC1_gain_table],
            "gainfield": [str(target_field)],
            "interp": cfg.casa_sc1_applycal.interp,
            "spwmap": cfg.casa_sc1_applycal.spwmap,
            "parang": cfg.casa_sc1_applycal.parang,
        }
        # Diagnostic: Gain amplitude vs time
        plot_gamp_time = {
            "msname": cfg.obs.msfile,
            "xaxis": cfg.casa_sc1_gainamp_time.xaxis,
            "yaxis": cfg.casa_sc1_gainamp_time.yaxis,
            "field": str(target_field),
            "iteraxis": cfg.casa_sc1_gainamp_time.iteraxis,
            "avgchannel": cfg.casa_sc1_gainamp_time.avgchannel,
            "avgtime": cfg.casa_sc1_gainamp_time.avgtime,
            "coloraxis": cfg.casa_sc1_gainamp_time.coloraxis,
            "expformat": cfg.casa_sc1_gainamp_time.expformat,
            "exprange": cfg.casa_sc1_gainamp_time.exprange,
            "plotfile": cfg.obs.basename +
                        plot_name[target.name] + "_" +
                        "sc1_gainamp_time.png",
            "overwrite": True,
        }
        # Diagnostic: Gain phase vs time
        plot_gphase_time = {
            "msname": cfg.obs.msfile,
            "xaxis": cfg.casa_sc1_gainphase_time.xaxis,
            "yaxis": cfg.casa_sc1_gainphase_time.yaxis,
            "field": str(target_field),
            "iteraxis": cfg.casa_sc1_gainphase_time.iteraxis,
            "avgchannel": cfg.casa_sc1_gainphase_time.avgchannel,
            "avgtime": cfg.casa_sc1_gainphase_time.avgtime,
            "coloraxis": cfg.casa_sc1_gainphase_time.coloraxis,
            "expformat": cfg.casa_sc1_gainphase_time.expformat,
            "exprange": cfg.casa_sc1_gainphase_time.exprange,
            "plotfile": cfg.obs.basename +
                        plot_name[target.name] + "_" +
                        "sc1_gainphase_time.png",
            "overwrite": True,
        }
        # Diagnostic: Gain amplitude vs freq
        plot_gamp_freq = {
            "msname": cfg.obs.msfile,
            "xaxis": cfg.casa_sc1_gainamp_freq.xaxis,
            "yaxis": cfg.casa_sc1_gainamp_freq.yaxis,
            "field": str(target_field),
            "iteraxis": cfg.casa_sc1_gainamp_freq.iteraxis,
            "avgchannel": cfg.casa_sc1_gainamp_freq.avgchannel,
            "avgtime": cfg.casa_sc1_gainamp_freq.avgtime,
            "coloraxis": cfg.casa_sc1_gainamp_freq.coloraxis,
            "expformat": cfg.casa_sc1_gainamp_freq.expformat,
            "exprange": cfg.casa_sc1_gainamp_freq.exprange,
            "plotfile": cfg.obs.basename +
                        plot_name[target.name] + "_" +
                        "sc1_gainamp_freq.png",
            "overwrite": True,
        }
        # Diagnostic: Gain phase vs freq
        plot_gphase_freq = {
            "msname": cfg.obs.msfile,
            "xaxis": cfg.casa_sc1_gainphase_freq.xaxis,
            "yaxis": cfg.casa_sc1_gainphase_freq.yaxis,
            "field": str(target_field),
            "iteraxis": cfg.casa_sc1_gainphase_freq.iteraxis,
            "avgchannel": cfg.casa_sc1_gainphase_freq.avgchannel,
            "avgtime": cfg.casa_sc1_gainphase_freq.avgtime,
            "coloraxis": cfg.casa_sc1_gainphase_freq.coloraxis,
            "expformat": cfg.casa_sc1_gainphase_freq.expformat,
            "exprange": cfg.casa_sc1_gainphase_freq.exprange,
            "plotfile": cfg.obs.basename +
                        plot_name[target.name] + "_" +
                        "sc1_gainphase_freq.png",
            "overwrite": True,
        }
        recipe.add("cab/casa_gaincal", "sc1_gain_calibration",
                   gain_cal_opts,
                   input=INPUT, output=OUTPUT,
                   label="casa_sc1_gaincal_%d:: Gain calibration" % target_field)

        recipe.add("cab/casa_applycal", "sc1_apply_calibration",
                   apply_cal_opts,
                   input=INPUT, output=OUTPUT,
                   label="casa_sc1_applycal_%d:: Apply calibration solutions to target" % target_field)
        recipe.add("cab/casa_plotms", "plot_gamp_time_sc1",
                   plot_gamp_time,
                   input=INPUT, output=OUTPUT,
                   label="plot_gamp_time_sc1_%d:: Gain amps with time" % target_field)
        recipe.add("cab/casa_plotms", "plot_gphase_time_sc1",
                   plot_gphase_time,
                   input=INPUT, output=OUTPUT,
                   label="plot_gphase_time_sc1_%d:: Gain phase with time" % target_field)
        recipe.add("cab/casa_plotms", "plot_gamp_freq_sc1",
                   plot_gamp_freq,
                   input=INPUT, output=OUTPUT,
                   label="plot_gamp_freq_sc1_%d:: Gain amps with freq" % target_field)
        recipe.add("cab/casa_plotms", "plot_gphase_freq_sc1",
                   plot_gphase_freq,
                   input=INPUT, output=OUTPUT,
                   label="plot_gphase_freq_sc1_%d:: Gain phase with freq" % target_field)
    recipe.run(["casa_sc1_gaincal_%d" % tf for tf in target_fields] +
               ["casa_sc1_applycal_%d" % tf for tf in target_fields] +
               ["plot_gamp_time_sc1_%d" % tf for tf in target_fields] +
               ["plot_gphase_time_sc1_%d" % tf for tf in target_fields] +
               ["plot_gamp_freq_sc1_%d" % tf for tf in target_fields] +
               ["plot_gphase_freq_sc1_%d" % tf for tf in target_fields])
