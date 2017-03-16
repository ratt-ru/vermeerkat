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
import math

def launch(cfg, INPUT, MSDIR, OUTPUT, **kwargs):
    plot_name = kwargs["plot_name"]
    targets = kwargs["targets"]
    target_fields = kwargs["target_fields"]
    recipe = stimela.Recipe("Phase only selfcal imaging engine", ms_dir=MSDIR)
    for target_field, target in zip(target_fields, targets):
        # make another mfs image and clean down gradually using masks
        model_prefix = cfg.obs.basename + "_LSM0_" + plot_name[target.name]
        model_name = model_prefix + ".lsm.html"
        imname_prefix = cfg.obs.basename + "_SC0_0_" + plot_name[target.name]
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC0_0_%d" % target_field,
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_sc0_0.column,
                       "weight": "briggs %2.f" % (cfg.wsclean_image_sc0_0.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_sc0_0.clean_iterations,
                       "mgain": cfg.wsclean_image_sc0_0.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_sc0_0.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_sc0_0.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_sc0_0.autothreshold,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_SC0_0_%d::wsclean" % target_field)
        maskname = imname_prefix + "_0_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask_sc0_0_%d" % target_field,
                   {
                       "image": '%s:output' % imname,
                       "output": maskname,
                       "sigma": cfg.cleanmask_sc0_0.sigma,
                       "iters": cfg.cleanmask_sc0_0.iters,
                       "boxes": cfg.cleanmask_sc0_0.kernel,
                       "no-negative": cfg.cleanmask_sc0_0.no_negative,
                       "tolerance" : cfg.cleanmask_sc0_0.tolerance,
                       "overlap" : cfg.cleanmask_sc0_0.overlap,
                   },
                   input=INPUT, output=OUTPUT,
                   label="make_clean_mask_sc0_0_%d::Make clean mask" % target_field)

        imname_prefix = cfg.obs.basename + "_SC0_1_" + plot_name[target.name]
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC0_1_%d" % target_field,
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_sc0_1.column,
                       "weight": "briggs %2.f" % (cfg.wsclean_image_sc0_1.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_sc0_1.clean_iterations,
                       "mgain": cfg.wsclean_image_sc0_1.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_sc0_1.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_sc0_1.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_sc0_1.autothreshold,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_SC0_1_%d::wsclean" % target_field)
        maskname = imname_prefix + "_1_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask_sc0_1_%d" % target_field,
                   {
                       "image": '%s:output' % imname,
                       "output": maskname,
                       "sigma": cfg.cleanmask_sc0_1.sigma,
                       "iters": cfg.cleanmask_sc0_1.iters,
                       "boxes": cfg.cleanmask_sc0_1.kernel,
                       "no-negative": cfg.cleanmask_sc0_1.no_negative,
                       "tolerance" : cfg.cleanmask_sc0_1.tolerance,
                       "overlap" : cfg.cleanmask_sc0_1.overlap,
                   },
                   input=INPUT, output=OUTPUT,
                   label="make_clean_mask_sc0_1_%d::Make clean mask" % target_field)

        imname_prefix = cfg.obs.basename + "_SC0_2_" + plot_name[target.name]
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC0_2_%d" % target_field,
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_sc0_2.column,
                       "weight": "briggs %2.f" % (cfg.wsclean_image_sc0_1.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_sc0_2.clean_iterations,
                       "mgain": cfg.wsclean_image_sc0_2.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_sc0_2.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_sc0_2.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_sc0_2.autothreshold,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_SC0_2_%d::wsclean" % target_field)
        maskname = imname_prefix + "_2_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask_sc0_2_%d" % target_field,
                   {
                       "image": '%s:output' % imname,
                       "output": maskname,
                       "sigma": cfg.cleanmask_sc0_2.sigma,
                       "iters": cfg.cleanmask_sc0_2.iters,
                       "boxes": cfg.cleanmask_sc0_2.kernel,
                       "no-negative": cfg.cleanmask_sc0_2.no_negative,
                       "tolerance" : cfg.cleanmask_sc0_2.tolerance,
                       "overlap" : cfg.cleanmask_sc0_2.overlap,
                   },
                   input=INPUT, output=OUTPUT,
                   label="make_clean_mask_sc0_2_%d::Make clean mask" % target_field)

        imname_prefix = cfg.obs.basename + "_SC0_3_" + plot_name[target.name]
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC0_3_%d" % target_field,
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_sc0_3.column,
                       "weight": "briggs %2.f" % (cfg.wsclean_image_sc0_3.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_sc0_3.clean_iterations,
                       "mgain": cfg.wsclean_image_sc0_3.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_sc0_3.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_sc0_3.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_sc0_3.autothreshold,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_SC0_3_%d::wsclean" % target_field)

    recipe.run(["image_SC0_0_%d" % t for t in target_fields] +
               ["make_clean_mask_sc0_0_%d" % t for t in target_fields] +
               ["image_SC0_1_%d" % t for t in target_fields] +
               ["make_clean_mask_sc0_1_%d" % t for t in target_fields] +
               ["image_SC0_2_%d" % t for t in target_fields] +
               ["make_clean_mask_sc0_2_%d" % t for t in target_fields] +
               ["image_SC0_3_%d" % t for t in target_fields])
