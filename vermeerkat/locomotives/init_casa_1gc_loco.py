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

import vermeerkat
import vermeerkat.dispatch_crew.caltables as vmct


def launch(cfg, INPUT, MSDIR, OUTPUT, **kwargs):
    bandpass_cal = kwargs["bandpass_cal"]
    calibrator_db = kwargs["calibrator_db"]
    bpcal_field = kwargs["bpcal_field"]
    gaincal_field = kwargs["gaincal_field"]
    bpcal_sol_int = kwargs["bpcal_sol_int"]
    target_fields = kwargs["target_fields"]

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
            "msname": cfg.obs.msfile,
            "field": str(bpcal_field),
            "standard": cfg.setjy_manual.standard,
            "fluxdensity": I,
            "spix": [a, b, c, d],
            "reffreq": "%.2fGHz" % (cfg.obs.freq_0 / 1e9),
            "usescratch": cfg.setjy_manual.usescratch,
            "scalebychan": cfg.setjy_manual.scalebychan,
            "spw": cfg.setjy_manual.spw,
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
            "msname": cfg.obs.msfile,
            "standard": cfg.setjy_auto.standard,
            "field": str(bpcal_field)
        }

    # @mauch points out some antenna positions may
    # be slightly off. This will cause a noticible
    # decorrelation on the longest baselines, so
    # better calibrate for this one (we can use the
    # bright bandpass and gain calibrator sources).
    delaycal_opts = {
        "msname": cfg.obs.msfile,
        "field": str(gaincal_field),
        "gaintype": cfg.delaycal.gaintype,
        "solint": cfg.delaycal.solint,
        "minsnr": cfg.delaycal.minsnr,
        "refant": cfg.obs.refant,
        "caltable": cfg.obs.delaycal_table,
        "calmode": cfg.delaycal.calmode,
    }

    # The bandpass calibrator may vary in phase over time
    # so lets do a preliminary phase cal to correct for this
    bp_phasecal_opts = {
        "msname": cfg.obs.msfile,
        "caltable": cfg.obs.phasecal_table,
        "field": str(bpcal_field),
        "refant": cfg.obs.refant,
        "calmode": cfg.phase0.calmode,
        "solint": cfg.obs.gain_sol_int,
        "solnorm": cfg.phase0.solnorm,
        "minsnr": cfg.phase0.minsnr,
        "gaintable": [cfg.obs.delaycal_table],
        "uvrange":cfg.phase0.uvrange,
    }

    # Then a bandpass calibration, lets work out a solution
    # per scan interval to see if things vary signifcantly
    # over time... this will be very bad for your reduction.
    bandpass_cal_opts = {
        "msname": cfg.obs.msfile,
        "caltable": cfg.obs.bpasscal_table,
        "field": str(bpcal_field),
        "spw": cfg.bandpass.spw,
        "refant": cfg.obs.refant,
        "solnorm": cfg.bandpass.solnorm,
        "combine": cfg.bandpass.combine,
        "solint": str(bpcal_sol_int) + "s",
        "bandtype": cfg.bandpass.bandtype,
        "minblperant": cfg.bandpass.minblperant,
        "minsnr": cfg.bandpass.minsnr,
        "gaintable": [cfg.obs.delaycal_table,
                      cfg.obs.phasecal_table],
        "interp": cfg.bandpass.interp,
        "uvrange": cfg.bandpass.uvrange,
    }

    # Finally we do a second order correction on the gain
    # cal source that is closest to the target fields
    gain_cal_opts = {
        "msname": cfg.obs.msfile,
        "caltable": cfg.obs.ampcal_table,
        "field": ",".join([str(x) for
                           x in [bpcal_field, gaincal_field]]),
        "spw": cfg.gaincal.spw,
        "solint": cfg.obs.gain_sol_int,
        "refant": cfg.obs.refant,
        "gaintype": cfg.gaincal.gaintype,
        "calmode": cfg.gaincal.calmode,
        "solnorm": cfg.gaincal.solnorm,
        "gaintable": [cfg.obs.delaycal_table,
                      cfg.obs.phasecal_table,
                      cfg.obs.bpasscal_table],
        "interp": cfg.gaincal.interp,
        "uvrange":cfg.gaincal.uvrange,
    }

    # Scale the gaincal solutions amplitude to that of the
    # bandpass calibrator (the flux scale reference)
    flux_scale_opts = {
        "msname": cfg.obs.msfile,
        "caltable": cfg.obs.ampcal_table,
        "fluxtable": cfg.obs.fluxcal_table,
        "reference": str(bpcal_field),
        "transfer": str(gaincal_field),
        "incremental": cfg.fluxscale.incremental,
    }

    # Apply gain solutions to all fields including
    # the calibrators so that we can diagnose problems more
    # easily. Later steps depend on this so don't remove
    # the gain application to calibrators
    apply_cal_opts = {
        "msname": cfg.obs.msfile,
        "field": ",".join([str(x) for x in
                           target_fields + [bpcal_field, gaincal_field]]),
        "gaintable": [cfg.obs.delaycal_table,
                      cfg.obs.phasecal_table,
                      cfg.obs.bpasscal_table,
                      cfg.obs.fluxcal_table],
        "gainfield": [str(x) for x in
                      gaincal_field,
                      bpcal_field,
                      bpcal_field,
                      gaincal_field],
        "interp": cfg.applycal.interp,
        "spwmap": cfg.applycal.spwmap,
        "parang": cfg.applycal.parang,
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

    # Go initial 1GC and further flagging
    init_1gc = ["setjy",
                "delaycal",
                "phase0",
                "bandpass",
                "gaincal",
                "fluxscale",
                "applycal"]
    recipe.run(init_1gc)
