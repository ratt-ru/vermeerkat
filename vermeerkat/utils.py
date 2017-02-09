from collections import namedtuple, defaultdict, OrderedDict
from pprint import pformat
import os

import numpy as np
import katdal

import vermeerkat
import vermeerkat.caltables as vmct

Scan = namedtuple("Scan", ["scan_index", "name", "tags", "radec", "length"])

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
    return [Scan(scan_index, target.name,
                        target.tags, target.radec(),
                        # timestamps changes depending on the scans
                        d.timestamps[-1] - d.timestamps[0])
             for (scan_index, state, target) in d.scans()
             if state == "track"]

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
    return [sum(s.length for s in field_scan_map[st.name])
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
    Select index and bandpass calibrator from a list
    """

    # See CASA setjy documentation for reasoning here

    # Base for Perley-Butler 2010 and 2013
    PERLEY_BUTLER_BASE = { "3C48",  "3C138", "3C147",
                            "3C196", "3C286", "3C295" }

    REYNOLDS_1994 = {"PKS 1934-638"}
    # PB2013 adds PKS 1934-638
    PERLEY_BUTLER_2010 = PERLEY_BUTLER_BASE.union(REYNOLDS_1994)
    # PB2013 adds 3C123
    PERLEY_BUTLER_2013 = PERLEY_BUTLER_BASE.union({"3C123"})
    SOUTHERN = set(vmct.calibrator_database().keys())

    # List of bandpass calibrators that we should prefer in order of priority.
    # Logic here is that there are only a few bandpass calibrators
    # See https://github.com/ska-sa/vermeerkat/issues/34
    prioritised_bpcals = [
        ("Perley-Butler 2013", PERLEY_BUTLER_2013),
        ("Perley-Butler 2010", PERLEY_BUTLER_2010),
        # Custom
        ("Southern", SOUTHERN),
    ]

    matches = None

    # Try for prioritised bandpass calibrators first
    for priority, (standard, candidates) in enumerate(prioritised_bpcals):
        matches = [(i, b, bp_scan_totals[i]) for i, b in enumerate(bpcals)
                                                    if b.name in candidates]

        if len(matches) > 0:
            index, bpcal, length = max(matches, key=lambda (i, bp, l): l)

            return index, bpcal, standard

    raise ValueError("Unable to find a bandpass calibrator in any standard. "
                    "Potential bandpass calibrators: '{}'. "
                    "Supported standards:\n{}".format(
                        bpcals, pformat(prioritised_bpcals)))