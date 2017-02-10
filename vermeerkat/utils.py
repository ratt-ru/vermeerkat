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
    bandpasses = []
    gaincals = []
    delaycals = []
    targets = []
    field_index = {}

    for scan in scans:
        categorise_field = False

        if scan.name in field_index:
            continue

        if "delaycal" in scan.tags:
            delaycals.append(scan)
            categorise_field = True

        if "bpcal" in scan.tags:
            bandpasses.append(scan)
            categorise_field = True

        if "gaincal" in scan.tags:
            gaincals.append(scan)
            categorise_field = True

        if "target" in scan.tags:
            targets.append(scan)
            categorise_field = True

        if not categorise_field:
            vermeerkat.log.warn("Not using observed field %s" % scan.name)
            continue

        field_index[scan.name] = len(field_index)

    return field_index, bandpasses, gaincals, delaycals, targets

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

def select_bandpass_calibrator(bpcals, bp_scan_totals):
    """
    Select index and bandpass calibrator from a list
    """

    if len(bpcals) == 0:
        raise ValueError("Zero bandpass calibrators supplied for selection")

    # See CASA setjy documentation for reasoning here

    # Base for Perley-Butler 2010 and 2013
    PERLEY_BUTLER_BASE = { "3C48",  "3C138", "3C147",
                            "3C196", "3C286", "3C295" }

    REYNOLDS_1994 = {"1934-638", "PKS 1934-638"}
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
        matches = [(i, b, st) for i, (b, st)
                    in enumerate(zip(bpcals, bp_scan_totals))
                    if b.name in candidates]

        if len(matches) > 0:
            break

    # Couldn't find anything, Southern should be last and will force
    # the main script to try and manually set the bandpass parameters
    # in the main script
    if matches is None:
        assert standard == "Southern"
        vermeerkat.log.info("Couldn't find any of the given calibrators '{}' "
                            "in our prioritised standards. "
                            "Selecting candidate with "
                            "longest observation time.".format(bpcals))

        matches = [(i, b, st) for i, (b, st)
                    in enumerate(zip(bpcals, bp_scan_totals))]

    # Choose bandpass with the longest observation time
    index, bpcal, length = max(matches, key=lambda (i, bp, l): l)
    return index, bpcal, standard


