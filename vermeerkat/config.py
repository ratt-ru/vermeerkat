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
import itertools
import os
import sys

import ruamel.yaml

import vermeerkat

__ARGS = None

def store_args(args):
    """ Store arguments for later retrieval """
    global __ARGS

    if __ARGS is not None:
        vermeerkat.log.warn("Replacing existing stored arguments '{}'"
                            "with '{}'.", __ARGS, args)

    __ARGS = args

def retrieve_args():
    """ Retrieve stored arguments """
    global __ARGS

    if __ARGS is None:
        raise ValueError("No arguments were stored. "
                        "Please call store_args first.")

    return __ARGS

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file '%s' does not exist!" % arg)

    return arg

def general_section_parser():
    """ Parses the general section """
    parser = argparse.ArgumentParser("General Section")

    # Solr server URL
    parser.add_argument('-s', '--solr-url',
        default='http://127.0.0.1/solr/core',
        help='URL of the Solr server')

    parser.add_argument('-f', '--hdf5-file',
        help='Name of the HDF5 file to download')

    parser.add_argument('-b', '--bandpass-calibrator',
        help="Imaging bandpass calibrator")

    parser.add_argument('-g', '--gain-calibrator',
        help="Imaging gain calibrator ")

    return parser

# Dictionary of argument parsers for particular sections
# Keys should correspond to associated sections in the
_ARGPARSERS = {
    'general': general_section_parser,
}

def configuration(args=None):
    """ Extract """

    #=========================================================
    # Handle the configuration file argument first,
    # if one is supplied use that for defaulting arguments
    # created further down the line, otherwise use the
    # default configuration file
    #=========================================================

    # Create parser object
    parser = argparse.ArgumentParser("VerMeerKAT")

    # Configuration file
    parser.add_argument('-c', '--config',
        type=lambda a: is_valid_file(parser, a),
        help='Configuration File',
        default=os.path.join(vermeerkat.install_path(),
            'conf', 'default.conf'),)

    # Parse configuration file arguments, if any
    args, remaining_args = parser.parse_known_args(args)

    # Lambda for transforming sections and options
    xformer = lambda s: s.lower().replace('-', '_')

    # Load in configuration options from the configuration file
    vermeerkat.log.info("Loading defaults from {}".format(args.config))

    with open(args.config, 'r') as f:
        file_config = ruamel.yaml.load(f, ruamel.yaml.SafeLoader)

    def _section_args(sections, args):
        """ Yields argparse.Namespace() containing arguments for each section """

        for section, (cfg_section, parser_factory) in sections.iteritems():
            # Extract configuration file options for this section
            # and use them as defaults for returning results

            section_defaults = file_config.get(cfg_section, None)

            # cfg_section might also be None...
            if section_defaults is None:
                section_defaults = {}

            # Transform keys
            section_defaults = { xformer(k): v for k, v
                                in section_defaults.iteritems() }

            # If present, use the argument parser to parse
            # (and validate) arguments
            if parser_factory:
                parser = parser_factory()
                parser.set_defaults(**section_defaults)
                section_args, args = parser.parse_known_args(args)
            # Otherwise just dump the section options into a namedtuple
            # that looks like one produced by argparse
            else:
                # vermeerkat.log.warn("No argument parser "
                #     "exists for Section '{}'. Options will "
                #     "not be validated for this section.".format(cfg_section))

                section_args = argparse.Namespace(**section_defaults)

            yield section, section_args

        if len(args) > 0:
            vermeerkat.log.warn("'{}' arguments were not parsed".format(args))

    # Find the list of sections that should be parsed.
    # Obtained from sections in the config file and sections in
    # the _ARGPARSERS dictionary.
    cfg_sections = { xformer(s): s for s in file_config.keys() }
    argp_sections = { xformer(k): pf for k, pf in _ARGPARSERS.iteritems() }
    all_sections = set(itertools.chain(cfg_sections.iterkeys(),
        argp_sections.iterkeys()))

    # Match configuration and arg parser sections up
    # after transforming their keys
    sections = { s: (cfg_sections.get(s, None), argp_sections.get(s, None))
        for s in all_sections }

    # Create the namespace using sections as arguments
    return argparse.Namespace(**{s: a for s, a in _section_args(
        sections, remaining_args)})
