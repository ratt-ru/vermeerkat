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
import copy

def launch(cfg, INPUT, MSDIR, OUTPUT, **kwargs):
    bandpass_cal = kwargs["bandpass_cal"]
    bpcal_field = kwargs["bpcal_field"]
    gaincal_field = kwargs["gaincal_field"]
    plot_name = kwargs["plot_name"]
    gain_cal = kwargs["gain_cal"]
    recipe = stimela.Recipe("1GC Diagnostics Engine", ms_dir=MSDIR)

    # Diagnostic: amplitude vs uv dist of the bandpass calibrator
    # @mauch points out we expect all baselines to observe the same amplitude
    # for the point source-like bandpass calibrator
    plot_ampuvdist_opts = {
        "msname": cfg.obs.msfile,
        "xaxis": cfg.plot_ampuvdist.xaxis,
        "yaxis": cfg.plot_ampuvdist.yaxis,
        "xdatacolumn": cfg.plot_ampuvdist.xdatacolumn,
        "ydatacolumn": cfg.plot_ampuvdist.ydatacolumn,
        "field": str(bpcal_field),
        "iteraxis": cfg.plot_ampuvdist.iteraxis,
        "avgchannel": cfg.plot_ampuvdist.avgchannel,
        "avgtime": cfg.plot_ampuvdist.avgtime,
        "coloraxis": cfg.plot_ampuvdist.coloraxis,
        "expformat": cfg.plot_ampuvdist.expformat,
        "exprange": cfg.plot_ampuvdist.exprange,
        "plotfile": cfg.obs.basename + "_bp_" +
                    plot_name[bandpass_cal.name] + "_" +
                    "ampuvdist.png",
        "overwrite": True,
    }

    # Diagnostic: phase vs uv dist of the bandpass calibrator
    plot_phaseuvdist_opts = {
        "msname": cfg.obs.msfile,
        "xaxis": cfg.plot_phaseuvdist.xaxis,
        "yaxis": cfg.plot_phaseuvdist.yaxis,
        "xdatacolumn": cfg.plot_phaseuvdist.xdatacolumn,
        "ydatacolumn": cfg.plot_phaseuvdist.ydatacolumn,
        "field": str(bpcal_field),
        "iteraxis": cfg.plot_phaseuvdist.iteraxis,
        "avgchannel": cfg.plot_phaseuvdist.avgchannel,
        "avgtime": cfg.plot_phaseuvdist.avgtime,
        "coloraxis": cfg.plot_phaseuvdist.coloraxis,
        "expformat": cfg.plot_phaseuvdist.expformat,
        "exprange": cfg.plot_phaseuvdist.exprange,
        "plotfile": cfg.obs.basename + "_bp_" +
                    plot_name[bandpass_cal.name] + "_" +
                    "phaseuvdist.png",
        "overwrite": True,
    }

    # Diagnostic: amplitude vs phase of bp calibrator per antenna
    plot_phaseball_opts = {
        "msname": cfg.obs.msfile,
        "xaxis": cfg.plot_phaseball.xaxis,
        "yaxis": cfg.plot_phaseball.yaxis,
        "xdatacolumn": cfg.plot_phaseball.xdatacolumn,
        "ydatacolumn": cfg.plot_phaseball.ydatacolumn,
        "field": str(bpcal_field),
        "iteraxis": cfg.plot_phaseball.iteraxis,
        "avgchannel": cfg.plot_phaseball.avgchannel,
        "avgtime": cfg.plot_phaseball.avgtime,
        "coloraxis": cfg.plot_phaseball.coloraxis,
        "expformat": cfg.plot_phaseball.expformat,
        "exprange": cfg.plot_phaseball.exprange,
        "plotfile": cfg.obs.basename + "_bp_" +
                    plot_name[bandpass_cal.name] + "_" +
                    "phaseball.png",
        "overwrite": True,
    }

    # Diagnostic: amplitude vs frequency of bp calibrator
    plot_amp_freq_opts = {
        "msname": cfg.obs.msfile,
        "xaxis": cfg.plot_amp_freq.xaxis,
        "yaxis": cfg.plot_amp_freq.yaxis,
        "xdatacolumn": cfg.plot_amp_freq.xdatacolumn,
        "ydatacolumn": cfg.plot_amp_freq.ydatacolumn,
        "field": str(bpcal_field),
        "iteraxis": cfg.plot_amp_freq.iteraxis,
        "avgchannel": cfg.plot_amp_freq.avgchannel,
        "avgtime": cfg.plot_amp_freq.avgtime,
        "coloraxis": cfg.plot_amp_freq.coloraxis,
        "expformat": cfg.plot_amp_freq.expformat,
        "exprange": cfg.plot_amp_freq.exprange,
        "plotfile": cfg.obs.basename + "_bp_" +
                    plot_name[bandpass_cal.name] + "_" +
                    "band.png",
        "overwrite": True,
    }

    # Diagnostic: phase vs time of bp calibrator
    # @oms points out slopes in this will indicate problems with
    # digitizer reference timing / probably also any uncorrected
    # antenna positions
    plot_phase_time_opts = {
        "msname": cfg.obs.msfile,
        "xaxis": cfg.plot_phase_time.xaxis,
        "yaxis": cfg.plot_phase_time.yaxis,
        "xdatacolumn": cfg.plot_phase_time.xdatacolumn,
        "ydatacolumn": cfg.plot_phase_time.ydatacolumn,
        "field": str(bpcal_field),
        "iteraxis": cfg.plot_phase_time.iteraxis,
        "avgchannel": cfg.plot_phase_time.avgchannel,
        "avgtime": cfg.plot_phase_time.avgtime,
        "coloraxis": cfg.plot_phase_time.coloraxis,
        "expformat": cfg.plot_phase_time.expformat,
        "exprange": cfg.plot_phase_time.exprange,
        "plotfile": cfg.obs.basename + "_bp_" +
                    plot_name[bandpass_cal.name] + "_" +
                    "phasevtime.png",
        "overwrite": True,
    }
    # Diagnostic: phase vs time of bp calibrator
    # For similar purposes as phase vs freq
    plot_phase_freq_opts = {
        "msname": cfg.obs.msfile,
        "xaxis": cfg.plot_phase_freq.xaxis,
        "yaxis": cfg.plot_phase_freq.yaxis,
        "xdatacolumn": cfg.plot_phase_freq.xdatacolumn,
        "ydatacolumn": cfg.plot_phase_freq.ydatacolumn,
        "field": str(bpcal_field),
        "iteraxis": cfg.plot_phase_freq.iteraxis,
        "avgchannel": cfg.plot_phase_freq.avgchannel,
        "avgtime": cfg.plot_phase_freq.avgtime,
        "coloraxis": cfg.plot_phase_freq.coloraxis,
        "expformat": cfg.plot_phase_freq.expformat,
        "exprange": cfg.plot_phase_freq.exprange,
        "plotfile": cfg.obs.basename + "_bp_" +
                    plot_name[bandpass_cal.name] + "_" +
                    "phasevfreq.png",
        "overwrite": True,
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
    plot_ampuvdist_opts["plotfile"] = cfg.obs.basename + "_gc_" + \
                                      plot_name[gain_cal.name] + "_" + \
                                      "ampuvdist.png"
    recipe.add("cab/casa_plotms", "plot_amp_v_uv_dist_of_gc_calibrator",
               plot_ampuvdist_opts,
               input=INPUT, output=OUTPUT,
               label="plot_ampuvdist_gc:: Diagnostic plot of amplitude with uvdist")

    plot_phaseuvdist_opts["field"] = str(gaincal_field)
    plot_phaseuvdist_opts["plotfile"] = cfg.obs.basename + "_gc_" + \
                                        plot_name[gain_cal.name] + "_" + \
                                        "phaseuvdist.png"
    recipe.add("cab/casa_plotms", "plot_phase_v_uv_dist_of_gc_calibrator",
               plot_phaseuvdist_opts,
               input=INPUT, output=OUTPUT,
               label="plot_phaseuvdist_gc:: Diagnostic plot of phase with uvdist")

    plot_phaseball_opts["field"] = str(gaincal_field)
    plot_phaseball_opts["plotfile"] = cfg.obs.basename + "_gc_" + \
                                      plot_name[gain_cal.name] + "_" + \
                                      "phaseball.png"
    recipe.add("cab/casa_plotms", "plot_amp_v_phase_of_gc_calibrator",
               plot_phaseball_opts,
               input=INPUT, output=OUTPUT,
               label="plot_phaseball_gc:: Diagnostic plot of phaseball")

    plot_amp_freq_opts["field"] = str(gaincal_field)
    plot_amp_freq_opts["plotfile"] = cfg.obs.basename + "_gc_" + \
                                     plot_name[gain_cal.name] + "_" + \
                                     "band.png"
    recipe.add("cab/casa_plotms", "plot_amp_v_freq_of_gc_calibrator",
               plot_amp_freq_opts,
               input=INPUT, output=OUTPUT,
               label="plot_amp_freq_gc:: Diagnostic plot of band")

    plot_phase_time_opts["field"] = str(gaincal_field)
    plot_phase_time_opts["plotfile"] = cfg.obs.basename + "_gc_" + \
                                       plot_name[gain_cal.name] + "_" + \
                                       "phasevtime.png"

    plot_phase_time_opts["field"] = str(gaincal_field)
    plot_phase_time_opts["plotfile"] = cfg.obs.basename + "_gc_" + \
                                       plot_name[gain_cal.name] + "_" + \
                                       "phasevtime.png"
    recipe.add("cab/casa_plotms", "plot_phase_vs_time_of_gc_calibrator",
               plot_phase_time_opts,
               input=INPUT, output=OUTPUT,
               label="plot_phase_time_gc:: Diagnostic plot of phase with time")

    plot_phase_freq_opts["field"] = str(gaincal_field)
    plot_phase_freq_opts["plotfile"] = cfg.obs.basename + "_gc_" + \
                                       plot_name[gain_cal.name] + "_" + \
                                       "phasevfreq.png"
    recipe.add("cab/casa_plotms", "plot_phase_vs_freq_of_gc_calibrator",
               plot_phase_freq_opts,
               input=INPUT, output=OUTPUT,
               label="plot_phase_freq_gc:: Diagnostic plot of phase with freq")

    # Diagnostic only: image bandpass
    recipe.add("cab/wsclean", "wsclean_bandpass",
               {
                   "msname": cfg.obs.msfile,
                   "column": cfg.wsclean_bandpass.column,
                   "weight": "briggs %.2f" % (cfg.wsclean_bandpass.robust),
                   "npix": cfg.obs.im_npix,
                   "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                   "clean_iterations": cfg.wsclean_bandpass.clean_iterations,
                   "mgain": cfg.wsclean_bandpass.mgain,
                   "channelsout": cfg.obs.im_numchans,
                   "joinchannels": cfg.wsclean_bandpass.joinchannels,
                   "field": str(bpcal_field),
                   "name": cfg.obs.basename + "_bp_" + plot_name[bandpass_cal.name],
                   "pol": cfg.general.imaging_pol,
               },
               input=INPUT, output=OUTPUT,
               label="image_bandpass::wsclean")

    # Diagnostic only: image gaincal
    recipe.add("cab/wsclean", "wsclean_gain",
               {
                   "msname": cfg.obs.msfile,
                   "column": cfg.wsclean_gain.column,
                   "weight": "briggs %.2f" % (cfg.wsclean_gain.robust),
                   "npix": cfg.obs.im_npix,
                   "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                   "clean_iterations": cfg.wsclean_gain.clean_iterations,
                   "mgain": cfg.wsclean_gain.mgain,
                   "channelsout": cfg.obs.im_numchans,
                   "joinchannels": cfg.wsclean_gain.joinchannels,
                   "field": str(gaincal_field),
                   "name": cfg.obs.basename + "_gc_" + plot_name[gain_cal.name],
                   "pol": cfg.general.imaging_pol,
               },
               input=INPUT, output=OUTPUT,
               label="image_gain::wsclean")

    diagnostics_1gc = ["plot_ampuvdist_bp",
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
    diagnostics_1gc += ["image_bandpass", "image_gain"]  # diagnostic only
    recipe.run(diagnostics_1gc)
