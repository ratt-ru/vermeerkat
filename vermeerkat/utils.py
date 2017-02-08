from collections import namedtuple, defaultdict

import numpy as np
import katdal

import vermeerkat

ScanInfo = namedtuple("ScanInfo", ["scan_index",
    "name", "tags", "radec", "length"])

class ObservationProperties(object):
    def __init__(self, **kwargs):
        for k, v in kwargs.iteritems():
            setattr(self, k, v)

def create_observation_properties(obs_meta, cfg):
    #Observation properties
    refant = str(obs_meta["RefAntenna"])
    correlator_integration_time = obs_meta["DumpPeriod"]
    gain_sol_int = str(correlator_integration_time * 3) + "s"
    freq_0 = obs_meta["CenterFrequency"]
    nchans = obs_meta["NumFreqChannels"]
    chan_bandwidth = obs_meta["ChannelWidth"]
    lambda_min = (299792458.0 / (freq_0 + (nchans // 2) * chan_bandwidth))
    telescope_max_baseline = 4.2e3 # TODO: need to automate calculation of this value
    angular_resolution = np.rad2deg(lambda_min /
                                    telescope_max_baseline *
                                    1.220) * 3600 #nyquest rate in arcsecs
    fov = cfg.general.fov * 3600 # 1 deg is good enough to cover FWHM of beam at L-Band
    sampling = cfg.general.sampling
    if sampling>1:
        raise ValueErro('PSF sampling is > 1. Please check your config file.')

    im_npix = int(fov / angular_resolution / sampling)
    bw_per_image_slice = 100.0e6
    im_numchans = int(np.ceil(obs_meta["ChannelWidth"] * nchans / bw_per_image_slice))

    return ObservationProperties(refant=refant,
        correlator_integration_time=correlator_integration_time,
        gain_sol_int=gain_sol_int,
        freq_0=freq_0,
        nchans=nchans,
        chan_bandwidth=chan_bandwidth,
        lambda_min=lambda_min,
        telescope_max_baseline=telescope_max_baseline,
        angular_resolution=angular_resolution,
        fov=fov,
        sampling=sampling,
        im_npix=im_npix,
        bw_per_image_slice=bw_per_image_slice,
        im_numchans=im_numchans)

def load_scans(h5filename):
    """ Load scan info from an h5 file """
    d = katdal.open(h5filename)
    return [ScanInfo(scan_index, target.name,
                        target.tags, target.radec(),
                        # timestamps changes depending on the scans
                        d.timestamps[-1] - d.timestamps[0])
             for (scan_index, state, target) in d.scans()
             if state == "track"]

def categorise_sources(scans):
    """ Categorise sources """
    bandpass_cal_candidates = []
    gain_cal_candidates = []
    targets = []
    source_index = {}

    for scan_info in scans:
        categorised_source = False

        if scan_info.name in source_index:
            continue

        if "bpcal" in scan_info.tags:
            bandpass_cal_candidates.append(scan_info)
            categorised_source = True

        if "gaincal" in scan_info.tags:
            gain_cal_candidates.append(scan_info)
            categorised_source = True

        if "target" in scan_info.tags:
            targets.append(scan_info)
            categorised_source = True

        if not categorised_source:
            vermeerkat.log.warn("Not using observed source %s" % name)
            continue

        source_index[scan_info.name] = len(source_index)

    return source_index, bandpass_cal_candidates, gain_cal_candidates, targets

def create_scan_map(scans):
    """ Map from scan target name to list of scans containing target """
    scan_map = defaultdict(list)

    for scan in scans:
        scan_map[scan.name].append(scan)

    return scan_map

def total_scan_times(scan_map, scan_targets):
    """
    Return a list of total_scan_times
    for each target in scan_targets
    """
    return [sum(s.length for s in scan_map[st.name])
                                    for st in scan_targets]

def select_gain_calibrator(targets, gaincals):
    """
    Select index and gain calibrator closest to the mean centre
    of the observation targets
    """

    # Compute mean target position
    mean_target = np.mean([t.radec for t in targets], axis=0)

    # Select gaincal candidate closest to mean target position
    sqrd_dist = (mean_target - [g.radec for g in gaincals])**2
    lmdistances = np.sum(sqrd_dist, axis=1)
    index = np.argmin(lmdistances)

    return index, gaincals[index]

def select_bandpass_calibrator(bpcals, bp_scan_totals):
    """
    Select index and bandpass calibrator with the longest total observation time
    """
    # Choose the bandpass calibrator with the longest observation time
    index, (bpcal, length) = max(enumerate(zip(bpcals, bp_scan_totals)),
        key=lambda (i, (bp, l)): l)

    return index, bpcal