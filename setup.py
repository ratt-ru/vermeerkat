#!/usr/bin/env python

import os

try:
      from setuptools import setup
except ImportError as e:
      from distutils.core import setup

requirements = [
 'stimela >= 1.0.0',
 'astropy>=2.0.11',
 'stimela>=1.3.1',
 'numpy>=1.13.1',
 'scipy>=0.19.1',
 'nbconvert>=5.3.1',
 'aplpy>=1.1.1',
 'matplotlib>=2.1.0',
 'jupyter>=1.0.0',
 'python-casacore>=2.2.1',
 'curses-menu>=0.5.0',
 'lmfit>=0.9.8',
 'GPy>=1.9.2',
 'npyscreen @ git+https://github.com/bennahugo/npyscreen.git@6bff9d3b65879dc',
]

PACKAGE_NAME = 'vermeerkat'
__version__ = '2.1.0'

setup(name = PACKAGE_NAME,
      version = __version__,
      description = "MeerKAT VermeerKAT pipeline",
      author = "B. Hugo",
      author_email = "bhugo@ska.ac.za",
      url = "https://github.com/bennahugo/vermeerkat",
      packages=[PACKAGE_NAME],
      install_requires = requirements,
      include_package_data = True,
      python_requires=">=3.6",
      entry_points={
          'console_scripts': ['vermeerkat=vermeerkat.bin.vermeerkat:main'],
      },
      license="GNU GPL v2",
      classifiers=[
                  "Development Status :: 3 - Alpha",
                  "Intended Audience :: Science/Research",
                  "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
                  "Operating System :: POSIX :: Linux",
                  "Programming Language :: Python",
                  "Topic :: Scientific/Engineering :: Astronomy"
              ]
 )
