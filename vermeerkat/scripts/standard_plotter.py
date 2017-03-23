#!/usr/bin/python
import argparse
import os
import sys

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

import vermeerkat
from vermeerkat.dispatch_crew.caltable_parser import read_caltable

ref_table = os.path.dirname(
        os.path.abspath(vermeerkat.__file__)) + "/southern_calibrators.txt"
calibrator_db = read_caltable(ref_table)

parser = argparse.ArgumentParser(description="Plot and compare measured fluxes "
                                 "against standard version")
parser.add_argument("calibrator_name",
                    metavar="calibrator_name",
                    type=str,
                    help="Name of the calibrator source, as specified in the "
                    "standard")
parser.add_argument("image_name",
                    metavar="image_name",
                    type=str,
                    nargs="+",
                    help="One or more FITS images of the calibrator source")

args = parser.parse_args(sys.argv[1:])

if args.calibrator_name not in calibrator_db:
    raise RuntimeError("Calibrator source not in standard")

ag = calibrator_db[args.calibrator_name]["a_ghz"]
bg = calibrator_db[args.calibrator_name]["b_ghz"]
cg = calibrator_db[args.calibrator_name]["c_ghz"]
dg = calibrator_db[args.calibrator_name]["d_ghz"]
std_freqs = np.linspace(0.1,5,1024)
log_std_freqs = np.log10(std_freqs)
fig = plt.figure(figsize=(15, 7), dpi=80)
ax1 = fig.add_subplot(111)
plt.title(args.calibrator_name, y=1.08)
ax1.set_xlabel("log(frequency [Hz])")
ax1.set_ylabel("log(Jy)")
log_std_freq_vals = ag + \
                    bg * np.log10(std_freqs) + \
                    cg * np.log10(std_freqs) ** 2 + \
                    dg * np.log10(std_freqs) ** 3

ax1.plot(log_std_freqs + np.log10(1e9),
         log_std_freq_vals, "b")

for img in args.image_name:
    with fits.open(img) as f:
        center_freq = f[0].header["CRVAL3"]
        freq_step = f[0].header["CDELT3"]
        nfreq = f[0].header["NAXIS3"]
        center_index = f[0].header["CRPIX3"] - 1 #FORTRAN indexing
        freq = np.arange(-center_index, nfreq, freq_step) + center_freq
        for fi in xrange(nfreq):
            freq_val = np.max(f[0].data[0, fi, :, :])
            rms = np.sqrt(np.mean(f[0].data[0, fi, :100, :100] ** 2))
            ax1.plot(np.log10(freq / 1e9) + np.log10(1e9),
                     np.log10(freq_val),"rx")
            ax1.errorbar(np.log10(freq / 1e9) + np.log10(1e9),
                         np.log10(freq_val),
                         yerr=np.log10((rms+freq_val)/freq_val),
                         c="g")

new_tick_locations = np.linspace(ax1.get_xlim()[0],
                                 ax1.get_xlim()[1], 8)
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(map(lambda x: "%.2f" % x,
                        10**(new_tick_locations - np.log10(1e9)) * 1e3))
ax2.set_xlabel(r"Frequency (MHz)")

plt.show()


