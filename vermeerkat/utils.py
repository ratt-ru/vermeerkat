import datetime
from collections import defaultdict
import os

import numpy as np
import katdal
from katpoint.target import construct_radec_target

import vermeerkat

class ObservationProperties(object):
    def __init__(self, **kwargs):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)

def merge_observation_metadata(cfg, obs_metadata):
    """ Merge observation metadata into the obs field of the configuration """
    cfg.obs = obs = ObservationProperties()

    # Configure filenames
    obs.h5file = obs_metadata['Filename']
    obs.basename = os.path.splitext(obs.h5file)[0]
    obs.msfile = ''.join((obs.basename, '.ms'))

    # Calibrator tables
    obs.delaycal_table = "%s.K0:output" % obs.basename
    obs.phasecal_table = "%s.G0:output" % obs.basename
    obs.ampcal_table = "%s.G1:output" % obs.basename
    obs.fluxcal_table = "%s.fluxscale:output" % obs.basename
    obs.bpasscal_table = "%s.B0:output" % obs.basename

    #Observation properties
    obs.refant = str(obs_metadata["RefAntenna"])
    obs.correlator_integration_time = obs_metadata["DumpPeriod"]
    obs.gain_sol_int = str(obs.correlator_integration_time * 3) + "s"
    obs.freq_0 = obs_metadata["CenterFrequency"]
    obs.nchans = obs_metadata["NumFreqChannels"]
    obs.chan_bandwidth = obs_metadata["ChannelWidth"]
    obs.lambda_min = (299792458.0 / (obs.freq_0 + (obs.nchans // 2) * obs.chan_bandwidth))
    obs.telescope_max_baseline = 4.2e3 # TODO: need to automate calculation of this value
    obs.angular_resolution = np.rad2deg(obs.lambda_min /
                                    obs.telescope_max_baseline *
                                    1.220) * 3600 #nyquest rate in arcsecs
    obs.fov = cfg.general.fov * 3600 # 1 deg is good enough to cover FWHM of beam at L-Band
    obs.sampling = cfg.general.sampling
    if obs.sampling>1:
        raise ValueErro('PSF sampling is > 1. Please check your config file.')

    obs.im_npix = int(obs.fov / obs.angular_resolution / obs.sampling)
    obs.bw_per_image_slice = 100.0e6
    obs.im_numchans = int(np.ceil(obs_metadata["ChannelWidth"] * obs.nchans / obs.bw_per_image_slice))

def load_scans(h5filename):
    """ Load scan info from an h5 file """
    d = katdal.open(h5filename)

    def _modify_target(target, scan_index, scan_length):
        target.scan_index = scan_index
        target.scan_length = scan_length
        return target

    return [_modify_target(target, scan_index, d.timestamps[-1] - d.timestamps[0])
                                    for (scan_index, state, target) in d.scans()
                                    if state == "track"]

def fmt_seconds(seconds, format=None):
    """ Formats seconds into %Hh%Mm%Ss format """

    if format is None:
        format = '%Hh%Mm%Ss'

    return datetime.datetime.utcfromtimestamp(seconds).strftime(format)

def categorise_fields(scans):
    """
    Categorise fields into targets, gain calibrators and bandpass calibrators
    """
    bandpass_cal_candidates = []
    gain_cal_candidates = []
    targets = []
    source_index = {}

    for scan in scans:
        categorise_field = False

        if scan.name in source_index:
            continue

        if "bpcal" in scan.tags:
            bandpass_cal_candidates.append(scan)
            categorise_field = True

        if "gaincal" in scan.tags:
            gain_cal_candidates.append(scan)
            categorise_field = True

        if "target" in scan.tags:
            targets.append(scan)
            categorise_field = True

        if not categorise_field:
            vermeerkat.log.warn("Not using observed field %s" % name)
            continue

        source_index[scan.name] = len(source_index)

    return source_index, bandpass_cal_candidates, gain_cal_candidates, targets

def create_field_scan_map(scans):
    """ Map from scan field to list of scans containing field """
    field_scan_map = defaultdict(list)

    for scan in scans:
        field_scan_map[scan.name].append(scan)

    return field_scan_map

def total_scan_times(field_scan_map, scan_targets):
    """
    Return a list of total_scan_times
    for each target in scan_targets
    """
    return [sum(s.scan_length for s in field_scan_map[st.name])
                                    for st in scan_targets]

def select_gain_calibrator(cfg, targets, gaincals):
    """
    Select index and gain calibrator closest to the mean centre
    of the observation targets
    """

    default_gaincal = cfg.general.gain_calibrator

    if default_gaincal:
        # If a default gain calibrator has been specified,
        # filter out other gain calibrators.
        vermeerkat.log.info("Gain calibrator manually "
                            "set to '%s'." % default_gaincal)

        gaincal_names = [g.name for g in gaincals]

        try:
            index = gaincal_names.index(default_gaincal)
            return index, gaincals[index]
        except ValueError:
            raise ValueError("'%s' not in list of "
                            "gain calibrators '%s'" % (
                                default_gaincal, gaincal_names))

    # Compute mean target position.
    # TODO: This assumes all targets are close to each other
    mean_radec = np.mean([t.radec() for t in targets], axis=0)

    # TODO: Is this the correct antenna to use for separation
    # below?
    ant = targets[0].antenna

    # Create a fake katpoint target
    mean_target = construct_radec_target(*mean_radec)

    # Select gaincal candidate with least angular distance
    # to mean target position
    angular_distances = [mean_target.separation(g, antenna=ant)
                                                for g in gaincals]
    index = np.argmin(angular_distances)

    return index, gaincals[index]

def select_bandpass_calibrator(cfg, bpcals, bp_scan_totals):
    """
    Select index and bandpass calibrator with the longest total observation time
    """
    # Choose the bandpass calibrator with the longest observation time

    default_bpcal = cfg.general.bandpass_calibrator

    if default_bpcal:
        # If a default gain calibrator has been specified,
        # filter out other gain calibrators.
        vermeerkat.log.info("Bandpass calibrator manually "
                            "set to '%s'." % default_bpcal)

        bpcal_names = [g.name for g in bpcals]

        try:
            index = bpcal_names.index(default_bpcal)
            return index, bpcals[index]
        except ValueError:
            raise ValueError("'%s' not in list of "
                            "bandpass calibrators '%s'" % (
                                default_bpcal, bpcal_names))


    index, (bpcal, length) = max(enumerate(zip(bpcals, bp_scan_totals)),
        key=lambda (i, (bp, l)): l)

    return index, bpcal