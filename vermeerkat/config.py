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

import argparse
import ConfigParser
import os
import urlparse
import sys

import vermeerkat

def configuration(args=None):
    """ Extract """

    print sys.argv

    # Create parser object
    parser = argparse.ArgumentParser("VerMeerKAT")

    #=========================================================
    # Handle the configuration file argument first,
    # if one is supplied use that for defaulting arguments
    # created further down the line, otherwise use the
    # default configuration file
    #=========================================================

    def is_valid_file(parser, arg):
        if not os.path.exists(arg):
            parser.error("The file %s does not exist!" % arg)

        return arg

    # Configuration file
    parser.add_argument('-c', '--config',
        type=lambda a: is_valid_file(parser, a),
        help='Configuration File',
        default=os.path.join(vermeerkat.install_path(),
            'conf', 'default.conf'),)

    # Parse configuration file arguments, if any
    args, remaining_argv = parser.parse_known_args()

    # Load in configuration options from the configuration file
    vermeerkat.log.info("Loading defaults from {}".format(args.config))
    config_parser = ConfigParser.SafeConfigParser()
    # Convert dashes in options to underscores
    # as this is how argparse will treat them
    config_parser.optionxform = lambda o: o.replace('-', '_')
    # Read the configuration file and extract options from the General section
    config_parser.read(args.config)
    config_defaults = dict(config_parser.items('General'))

    #==========================================================
    # Handle the rest of the arguments
    # A default value should be supplied, even if it's
    # to just illustrate what the argument should look like
    #==========================================================

    # Solr server URL
    parser.add_argument('-s', '--solr-url',
        default='http://127.0.0.1/solr/core',
        help='URL of the Solr server')

    # Set any defaults taken from the configuration file
    parser.set_defaults(**config_defaults)
    # Parse the rest of the command line arguments
    return parser.parse_args(remaining_argv)
