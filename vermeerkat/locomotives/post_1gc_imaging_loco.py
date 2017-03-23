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
import stimela

def launch(cfg, INPUT, MSDIR, OUTPUT, **kwargs):
    targets = kwargs["targets"]
    field_index = kwargs["field_index"]
    target_fields = kwargs["target_fields"]
    recipe = stimela.Recipe("1GC Imaging Engine", ms_dir=MSDIR)
    # imaging
    for target_field, target in zip(target_fields, targets):
        imname_prefix = cfg.obs.basename + "_1GC_" +"0_" + \
            target.name.replace(" ","_")
        imname = imname_prefix + "-MFS-image.fits"
        # gradually clean down using masks
        recipe.add("cab/wsclean", "wsclean_image_1gc_0_%d" % field_index[target.name],
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_1gc_0.column,
                       "weight": "briggs %.2f" % (cfg.wsclean_image_1gc_0.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_1gc_0.clean_iterations,
                       "mgain": cfg.wsclean_image_1gc_0.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_1gc_0.joinchannels,
                       "field": target_field,
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_1gc_0.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_1gc_0.autothreshold,
                       "pol": cfg.general.imaging_pol,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_1gc_0_%d::wsclean" % target_field)
        maskname = imname_prefix + "_0_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask_1gc_0_%d" % target_field,
                   {
                       "image": '%s:output' % imname,
                       "output": maskname,
                       "sigma": cfg.cleanmask_1gc_0.sigma,
                       "iters": cfg.cleanmask_1gc_0.iters,
                       "boxes": cfg.cleanmask_1gc_0.kernel,
                       "no-negative": cfg.cleanmask_1gc_0.no_negative,
                       "tolerance" : cfg.cleanmask_1gc_0.tolerance,
                       "overlap" : cfg.cleanmask_1gc_0.overlap,
                   },
                   input=INPUT, output=OUTPUT,
                   label="make_clean_mask_1gc_0_%d::Make clean mask" % target_field)

        imname_prefix = cfg.obs.basename + "_1GC_" + "1_" + \
            target.name.replace(" ","_")
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_image_1gc_1_%d" % field_index[target.name],
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_1gc_1.column,
                       "weight": "briggs %.2f" % (cfg.wsclean_image_1gc_1.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_1gc_1.clean_iterations,
                       "mgain": cfg.wsclean_image_1gc_1.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_1gc_1.joinchannels,
                       "field": target_field,
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_1gc_1.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_1gc_1.autothreshold,
                       "fitsmask": "%s:output" % maskname,
                       "pol": cfg.general.imaging_pol,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_1gc_1_%d::wsclean" % target_field)

        maskname = imname_prefix + "_1_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask_1gc_1_%d" % target_field,
           {
               "image": '%s:output' % imname,
               "output": maskname,
               "sigma": cfg.cleanmask_1gc_1.sigma,
               "iters": cfg.cleanmask_1gc_1.iters,
               "boxes": cfg.cleanmask_1gc_1.kernel,
               "no-negative": cfg.cleanmask_1gc_1.no_negative,
               "tolerance" : cfg.cleanmask_1gc_1.tolerance,
               "overlap" : cfg.cleanmask_1gc_1.overlap,

           },
           input=INPUT, output=OUTPUT,
           label="make_clean_mask_1gc_1_%d::Make clean mask" % target_field)

        imname_prefix = cfg.obs.basename + "_1GC_" + "2_" + \
            target.name.replace(" ","_")
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_image_1gc_2_%d" % field_index[target.name],
               {
                   "msname": cfg.obs.msfile,
                   "column": cfg.wsclean_image_1gc_2.column,
                   "weight": "briggs %.2f" % (cfg.wsclean_image_1gc_2.robust),
                   "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                   "trim" : cfg.obs.im_npix,
                   "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                   "clean_iterations": cfg.wsclean_image_1gc_2.clean_iterations,
                   "mgain": cfg.wsclean_image_1gc_2.mgain,
                   "channelsout": cfg.obs.im_numchans,
                   "joinchannels": cfg.wsclean_image_1gc_2.joinchannels,
                   "field": target_field,
                   "name": imname_prefix,
                   "taper-gaussian": cfg.wsclean_image_1gc_2.taper_gaussian,
                   "auto-threshold": cfg.wsclean_image_1gc_2.autothreshold,
                   "fitsmask": "%s:output" % maskname,
                   "pol": cfg.general.imaging_pol,
               },
               input=INPUT, output=OUTPUT,
               label="image_1gc_2_%d::wsclean" % target_field)

        maskname = imname_prefix + "_2_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask_1gc_2_%d" % target_field,
               {
                   "image": '%s:output' % imname,
                   "output": maskname,
                   "sigma": cfg.cleanmask_1gc_2.sigma,
                   "iters": cfg.cleanmask_1gc_2.iters,
                   "boxes": cfg.cleanmask_1gc_2.kernel,
                   "no-negative": cfg.cleanmask_1gc_2.no_negative,
                   "tolerance" : cfg.cleanmask_1gc_2.tolerance,
                   "overlap" : cfg.cleanmask_1gc_2.overlap,
               },
               input=INPUT, output=OUTPUT,
               label="make_clean_mask_1gc_2_%d::Make clean mask" % target_field)

        imname_prefix = cfg.obs.basename + "_1GC_" + "3_" + \
            target.name.replace(" ","_")
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_image_1gc_3_%d" % field_index[target.name],
               {
                   "msname": cfg.obs.msfile,
                   "column": cfg.wsclean_image_1gc_3.column,
                   "weight": "briggs %.2f" % (cfg.wsclean_image_1gc_3.robust),
                   "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                   "trim" : cfg.obs.im_npix,
                   "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                   "clean_iterations": cfg.wsclean_image_1gc_3.clean_iterations,
                   "mgain": cfg.wsclean_image_1gc_3.mgain,
                   "channelsout": cfg.obs.im_numchans,
                   "joinchannels": cfg.wsclean_image_1gc_3.joinchannels,
                   "field": target_field,
                   "name": imname_prefix,
                   "taper-gaussian": cfg.wsclean_image_1gc_3.taper_gaussian,
                   "auto-threshold": cfg.wsclean_image_1gc_3.autothreshold,
                   "fitsmask": "%s:output" % maskname,
                   "pol": cfg.general.imaging_pol,
               },
               input=INPUT, output=OUTPUT,
               label="image_1gc_3_%d::wsclean" % target_field)

    recipe.run(["image_1gc_0_%d" % t for t in target_fields] +
               ["make_clean_mask_1gc_0_%d" % t for t in target_fields] +
               ["image_1gc_1_%d" % t for t in target_fields] +
               ["make_clean_mask_1gc_1_%d" % t for t in target_fields] +
               ["image_1gc_2_%d" % t for t in target_fields] +
               ["make_clean_mask_1gc_2_%d" % t for t in target_fields] +
               ["image_1gc_3_%d" % t for t in target_fields])
