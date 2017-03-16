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

import stimela

import vermeerkat
import vermeerkat.conf.config as vmc
import vermeerkat.dispatch_crew.caltables as vmct
import vermeerkat.dispatch_crew.observation as vmo
import vermeerkat.dispatch_crew.utils as vmu
from vermeerkat.locomotives import converter_loco, \
                                   flagging_loco, \
                                   init_casa_1gc_loco, \
                                   postcal_flagging_loco, \
                                   secondpass_casa_1gc_loco, \
                                   post_1gc_diagnostics_loco, \
                                   post_1gc_imaging_loco, \
                                   phase_selfcal_loco, \
                                   phaseamp_selfcal_loco, \
                                   post_p_selfcal_imaging_loco, \
                                   post_ap_selfcal_imaging_loco, \
                                   phase_selfcal_casa_loco, \
                                   phaseamp_selfcal_casa_loco

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
    if not cfg.general.skip_conversion:
        converter_loco.launch(cfg, INPUT, MSDIR, OUTPUT)

    #########################################################################
    #
    # Preliminary RFI, band and autocorr flagging
    #
    #########################################################################
    if not cfg.general.skip_rfi_flagging:
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
    if not cfg.general.skip_initial_1gc:
        init_casa_1gc_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                                  bandpass_cal=bandpass_cal,
                                  calibrator_db=calibrator_db,
                                  bpcal_field=bpcal_field,
                                  gaincal_field=gaincal_field,
                                  bpcal_sol_int=bpcal_sol_int,
                                  target_fields=target_fields)
    if not cfg.general.skip_secondpass_flagging:
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
    if not cfg.general.skip_1gc_recalibration:
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
    if not cfg.general.skip_1gc_diagnostics:
        post_1gc_diagnostics_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                                         bandpass_cal=bandpass_cal,
                                         bpcal_field=bpcal_field,
                                         gaincal_field=gaincal_field,
                                         plot_name=plot_name,
                                         gain_cal=gain_cal)

    #########################################################################
    #
    # Post 1GC imaging
    #
    #########################################################################
    if not cfg.general.skip_1gc_imaging:
        post_1gc_imaging_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                                     targets=targets,
                                     field_index=field_index,
                                     target_fields=target_fields)

    #########################################################################
    #
    # Phase only self calibration (2nd gen)
    #
    #########################################################################
    sc0_loco = phase_selfcal_loco if cfg.general.replace_casa_with_mt_selfcal else phase_selfcal_casa_loco
    if not cfg.general.skip_phaseonly_selfcal:
        sc0_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                        plot_name=plot_name,
                        targets=targets,
                        target_fields=target_fields)
    if not cfg.general.skip_phaseonly_selfcal_imaging:
        post_p_selfcal_imaging_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                                           targets=targets,
                                           field_index=field_index,
                                           target_fields=target_fields,
                                           plot_name=plot_name)

    #########################################################################
    #
    # Phase Amplitude self calibration (2nd gen)
    #
    # After most of the significant phase error has been corrected
    # we should be able to clean deeper, so rerun this time with amplitude
    #########################################################################
    sc1_loco = phaseamp_selfcal_loco if cfg.general.replace_casa_with_mt_selfcal else phaseamp_selfcal_casa_loco
    if not cfg.general.skip_ampphase_selfcal:
        sc1_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                        plot_name=plot_name,
                        targets=targets,
                        target_fields=target_fields)
    if not cfg.general.skip_ampphase_selfcal_imaging:
        post_ap_selfcal_imaging_loco.launch(cfg, INPUT, MSDIR, OUTPUT,
                                            targets=targets,
                                            field_index=field_index,
                                            target_fields=target_fields,
                                            plot_name=plot_name)

