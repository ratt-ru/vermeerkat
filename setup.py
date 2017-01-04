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

import os
import shutil
import subprocess
import sys
from setuptools import setup, find_packages

PACKAGE_NAME = 'vermeerkat'

def get_version():
    # Versioning code here, based on
    # http://blogs.nopcode.org/brainstorm/2013/05/20/pragmatic-python-versioning-via-setuptools-and-git-tags/

    # Fetch version from git tags, and write to version.py.
    # Also, when git is not available (PyPi package), use stored version.py.
    version_py = os.path.join(PACKAGE_NAME, 'version.py')

    try:
        version_git = subprocess.check_output(['git', 'describe', '--tags']).rstrip()
    except:
        with open(version_py, 'r') as fh:
            version_git = open(version_py).read().strip().split('=')[-1].replace('"','')

    version_msg = "# Do not edit this file, pipeline versioning is governed by git tags"

    with open(version_py, 'w') as fh:
        fh.write(version_msg + os.linesep + "__version__=\"" + version_git +"\"")

    return version_git

def readme():
    with open('README.md') as f:
        return f.read()

setup(name=PACKAGE_NAME,
    version=get_version(),
    description='MeerKAT Imaging Pipeline',
    long_description=readme(),
    url='http://github.com/SpheMakh/VerMeerKAT',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    author='Simon Perkins',
    author_email='simon.perkins@gmail.com',
    license='GPL2',
    packages=find_packages(),
    install_requires=[
        'numpy >= 1.11.2',
        'progressbar2 >= 3.11.0',
        'pysolr >= 3.6.0',
        'requests >= 2.12.1',
        'stimela >= 0.2.0', # Force install of master branch
        'katdal >= 0.7',
        'numpy >= 1.11.2'
    ],
    dependency_links=[
        'https://github.com/SpheMakh/Stimela/tarball/master#egg=stimela-0.2.0'
    ],
    package_data={ PACKAGE_NAME: 'conf/*.conf' },
    include_package_data=True,
    zip_safe=False)
