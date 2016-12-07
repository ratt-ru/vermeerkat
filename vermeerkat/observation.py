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

import collections
import itertools
import os
import sys

import vermeerkat

def query_recent_observations(solr_url):
    """
    Find recent telescope observations suitable for imaging
    """
    import pysolr

    ONE_MINUTE = 60
    ONE_HOUR = 60 * ONE_MINUTE

    # Create a quick field+query type
    fq = collections.namedtuple('field_query', 'field query')

    query_list = [
        # Only want MeerKAT AR1 Telescope Products
        fq('CAS.ProductTypeName', 'MeerKATAR1TelescopeProduct'),
        # Only want 4K mode
        fq('NumFreqChannels', '4096'),
        # Observations from the last 3 days
        fq('StartTime', '[NOW-3DAYS TO NOW]')]

    # Construct the query
    search = ' AND '.join('%s:%s' % (fq.field, fq.query) for fq in query_list)

    vermeerkat.log.info("Querying solr server '%s' "
                        "with query string '%s'." % (solr_url, search))

    archive = pysolr.Solr(solr_url)

    res = archive.search(search, sort='StartTime desc', rows=1000)

    def _observation_filter(solr_result):
        """
        Filter out KAT observations that don't cut the mustard
        for imaging purposes.
        """

        # Create an observation string for logging
        filename, description = (solr_result[f] for f in
            ('Filename', 'Description'))
        observation = '{} {}'.format(filename, description[:50])

        # Its only worth imaging observations with
        # gain and bandpass calibrations
        targets = solr_result.get('KatpointTargets', [])
        gaincals = sum('gaincal' in s for s in targets)
        bandpasses = sum('bpcal' in s for s in targets)

        if gaincals == 0:
            vermeerkat.log.warn('Ignoring "{}", no gain calibrators '
                                'are present.'.format(observation))
            return False

        if bandpasses == 0:
            vermeerkat.log.warn('Ignoring "{}", no band pass calibrators '
                                ' are present.'.format(observation))
            return False

        # Don't take anything less than 2 hours in duration
        duration = solr_result['Duration']

        if duration <= 2*ONE_HOUR:
            vermeerkat.log.warn('Ignoring "{}", observation is '
                'less than 2 hours "{}s".'.format(observation, duration))
            return False

        return True

    # Filter observations
    return filter(_observation_filter, res)


def download_observation(observation, directory):
    """ Download the specified observation """
    import requests
    import progressbar

    ONE_KB = 1024
    ONE_MB = ONE_KB*ONE_KB

    targets = observation.get('KatpointTargets', [])
    gaincals = sum('gaincal' in s for s in targets)
    bandpasses = sum('bpcal' in s for s in targets)

    vermeerkat.log.info('%s -- %-2d -- %-2d-- %s -- %s -- %s' % (
        observation['Filename'], gaincals, bandpasses,
        observation['StartTime'], observation['Description'], observation['Duration']))

    # Infer the HTTP location from the KAT archive file location
    location = observation['FileLocation'][0]
    location = location.replace('/var/kat', 'http://kat-archive.kat.ac.za', 1)
    filename = observation['Filename']
    url = os.path.join(location, filename)

    filename = os.path.join(directory, filename)
    file_exists = os.path.exists(filename) and os.path.isfile(filename)
    local_file_size = os.path.getsize(filename) if file_exists else 0
    headers = { "Range" : "bytes={}-".format(local_file_size) }

    r = requests.get(url, headers=headers, stream=True)
    print r.headers

    # Server doesn't care about range requests and is just
    # sending the entire file
    if r.status_code == 200:
        vermeerkat.log.info("Downloading '{}'")
        remote_file_size = r.headers.get('Content-Length', None)
        file_exists = False
        local_file_size = 0
    elif r.status_code == 206:
        if local_file_size > 0:
            vermeerkat.log.info("'{}' already exists, "
                "resuming download from {}.".format(
                    filename, local_file_size))

        # Create a fake range if none exists
        fake_range = "{}-{}/{}".format(local_file_size, sys.maxint,
            sys.maxint - local_file_size)

        remote_file_size = r.headers.get('Content-Range', fake_range)
        remote_file_size = remote_file_size.split('/')[-1]
    elif r.status_code == 416:
        vermeerkat.log.info("'{}' already downloaded".format(filename))
        remote_file_size = local_file_size
        return filename
    else:
        raise ValueError("HTTP Error Code {}".format(r.status_code))

    vermeerkat.log.info('%s %s %s' % (url, remote_file_size, r.status_code))

    f = open(filename, 'ab' if file_exists else 'wb')
    bar = (progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        if remote_file_size is None else progressbar.ProgressBar(
            max_value=int(remote_file_size)))

    #Download chunks of file and write to disk
    try:
        with f, bar:
            downloaded = local_file_size
            for chunk in r.iter_content(chunk_size=ONE_MB):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    bar.update(downloaded)
    except KeyboardInterrupt as kbe:
        pass
        #log.warn("Quitting download on Keyboard Interrupt")

    return filename
