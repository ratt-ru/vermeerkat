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
import collections
import ConfigParser
import itertools
import os
import sys

import vermeerkat

def general_section_parser():
    """ Parses the general section """
    parser = argparse.ArgumentParser("General Section")

    # Solr server URL
    parser.add_argument('-s', '--solr-url',
        default='http://127.0.0.1/solr/core',
        help='URL of the Solr server')

    return parser

# Dictionary of argument parsers for particular sections
# Keys should correspond to associated sections in the
_ARGPARSERS = { 'general': general_section_parser }

def configuration(args=None):
    """ Extract """

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
    args, remaining_args = parser.parse_known_args(args)

    # Lambda for transforming sections and options
    xformer = lambda s: s.lower().replace('-', '_')

    # Load in configuration options from the configuration file
    vermeerkat.log.info("Loading defaults from {}".format(args.config))
    config_parser = ConfigParser.SafeConfigParser()
    # Convert dashes in options to underscores
    # as this is how argparse will treat them
    config_parser.optionxform = xformer
    # Read the configuration file and extract options from the General section
    config_parser.read(args.config)

    def _section_args(sections, args):
        """ Yields namedtuples containing arguments for each section """

        for section, (cfg_section, parser_factory) in sections.iteritems():
            # Extract configuration file options for this section
            # and use them as defaults for returning results

            section_defaults = ({} if cfg_section is None
                else dict(config_parser.items(cfg_section)))

            # If present, use the argument parser to parse
            # (and validate) arguments
            if parser_factory:
                parser = parser_factory()
                parser.set_defaults(**section_defaults)
                section_args, args = parser.parse_known_args(args)
            # Otherwise just dump the section options into a namedtuple
            # that looks like one produced by argparse
            else:
                vermeerkat.log.warn("No argument parser "
                    "exists for Section '{}'. Options will "
                    "not be validated for this section.".format(cfg_section))

                T = collections.namedtuple('Namespace',
                    section_defaults.keys())

                section_args = T(**section_defaults)

            yield section_args

    # Find the list of sections that should be parsed.
    # Obtained from sections in the config file and sections in
    # the _ARGPARSERS dictionary.
    cfg_sections = { xformer(s): s for s in config_parser.sections() }
    argp_sections = { xformer(k): pf for k, pf in _ARGPARSERS.iteritems() }
    all_sections = set(itertools.chain(cfg_sections.iterkeys(),
        argp_sections.iterkeys()))

    # Match configuration and arg parser sections up
    # after transforming their keys
    sections = { s: (cfg_sections.get(s, None), argp_sections.get(s, None))
        for s in all_sections }

    # Create a global namespace type that holds argparse
    # results for each section. Looks like an argparse args Namespace
    global_namespace = collections.namedtuple('Namespace', all_sections)

    # Create the namespace using sections as arguments
    return global_namespace(*(s for s in _section_args(
        sections, remaining_args)))