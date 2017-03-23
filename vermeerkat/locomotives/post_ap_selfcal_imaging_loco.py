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
import vermeerkat

def launch(cfg, INPUT, MSDIR, OUTPUT, **kwargs):
    plot_name = kwargs["plot_name"]
    targets = kwargs["targets"]
    target_fields = kwargs["target_fields"]
    do_full_restore = kwargs["do_full_restore"]
    recipe = stimela.Recipe("Phase amplitude selfcal imaging engine", ms_dir=MSDIR)
    for target_field, target in zip(target_fields, targets):
        # make another mfs image and clean down gradually using masks
        model_prefix = cfg.obs.basename + "_LSM1_" + plot_name[target.name]
        model_name = model_prefix + ".lsm.html"
        imname_prefix = cfg.obs.basename + "_SC1_0_" + plot_name[target.name]
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC1_0_%d" % target_field,
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_sc1_0.column,
                       "weight": "briggs %2.f" % (cfg.wsclean_image_sc1_0.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_sc1_0.clean_iterations,
                       "mgain": cfg.wsclean_image_sc1_0.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_sc1_0.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_sc1_0.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_sc1_0.autothreshold,
                       "pol": cfg.general.imaging_pol,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_SC1_0_%d::wsclean" % target_field)
        maskname = imname_prefix + "_0_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask_sc1_0_%d" % target_field,
                   {
                       "image": '%s:output' % imname,
                       "output": maskname,
                       "sigma": cfg.cleanmask_sc1_0.sigma,
                       "iters": cfg.cleanmask_sc1_0.iters,
                       "boxes": cfg.cleanmask_sc1_0.kernel,
                       "no-negative": cfg.cleanmask_sc1_0.no_negative,
                       "tolerance" : cfg.cleanmask_sc1_0.tolerance,
                       "overlap" : cfg.cleanmask_sc1_0.overlap,

                   },
                   input=INPUT, output=OUTPUT,
                   label="make_clean_mask_sc1_0_%d::Make clean mask" % target_field)

        imname_prefix = cfg.obs.basename + "_SC1_1_" + plot_name[target.name]
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC1_1_%d" % target_field,
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_sc1_1.column,
                       "weight": "briggs %2.f" % (cfg.wsclean_image_sc1_1.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_sc1_1.clean_iterations,
                       "mgain": cfg.wsclean_image_sc1_1.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_sc1_1.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_sc1_1.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_sc1_1.autothreshold,
                       "pol": cfg.general.imaging_pol,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_SC1_1_%d::wsclean" % target_field)
        maskname = imname_prefix + "_1_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask_sc1_1_%d" % target_field,
                   {
                       "image": '%s:output' % imname,
                       "output": maskname,
                       "sigma": cfg.cleanmask_sc1_1.sigma,
                       "iters": cfg.cleanmask_sc1_1.iters,
                       "boxes": cfg.cleanmask_sc1_1.kernel,
                       "no-negative": cfg.cleanmask_sc1_1.no_negative,
                       "tolerance" : cfg.cleanmask_sc1_1.tolerance,
                       "overlap" : cfg.cleanmask_sc1_1.overlap,
                   },
                   input=INPUT, output=OUTPUT,
                   label="make_clean_mask_sc1_1_%d::Make clean mask" % target_field)

        imname_prefix = cfg.obs.basename + "_SC1_2_" + plot_name[target.name]
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC1_2_%d" % target_field,
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_sc1_2.column,
                       "weight": "briggs %2.f" % (cfg.wsclean_image_sc1_1.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_sc1_2.clean_iterations,
                       "mgain": cfg.wsclean_image_sc1_2.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_sc1_2.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_sc1_2.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_sc1_2.autothreshold,
                       "pol": cfg.general.imaging_pol,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_SC1_2_%d::wsclean" % target_field)
        maskname = imname_prefix + "_2_MASK.fits"
        recipe.add("cab/cleanmask", "make_clean_mask_sc1_2_%d" % target_field,
                   {
                       "image": '%s:output' % imname,
                       "output": maskname,
                       "sigma": cfg.cleanmask_sc1_2.sigma,
                       "iters": cfg.cleanmask_sc1_2.iters,
                       "boxes": cfg.cleanmask_sc1_2.kernel,
                       "no-negative": cfg.cleanmask_sc1_2.no_negative,
                       "tolerance" : cfg.cleanmask_sc1_2.tolerance,
                       "overlap" : cfg.cleanmask_sc1_2.overlap,
                   },
                   input=INPUT, output=OUTPUT,
                   label="make_clean_mask_sc1_2_%d::Make clean mask" % target_field)

        imname_prefix = cfg.obs.basename + "_SC1_3_" + plot_name[target.name]
        imname = imname_prefix + "-MFS-image.fits"
        recipe.add("cab/wsclean", "wsclean_SC1_3_%d" % target_field,
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_image_sc1_3.column,
                       "weight": "briggs %2.f" % (cfg.wsclean_image_sc1_3.robust),
                       "npix": int(cfg.obs.im_npix * cfg.obs.padding),
                       "trim" : cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_image_sc1_3.clean_iterations,
                       "mgain": cfg.wsclean_image_sc1_3.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_image_sc1_3.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "taper-gaussian": cfg.wsclean_image_sc1_3.taper_gaussian,
                       "auto-threshold": cfg.wsclean_image_sc1_3.autothreshold,
                       "pol": cfg.general.imaging_pol,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_SC1_3_%d::wsclean" % target_field)

        # full restore into the final mfs
        if do_full_restore:
            imname_mfs_fullrest = imname_prefix + "-MFS-image.full_restore.fits"
            recipe.add("cab/tigger_restore", "full_restore_SC1_%d" % target_field,
                       {
                           "input-image": "%s:output" % imname,
                           "input-skymodel": "%s:output" % model_name,
                           "output-image": '%s:output' % imname_mfs_fullrest,
                           "freq": cfg.obs.freq_0,
                       },
                       input=INPUT, output=OUTPUT,
                       label="fullrestore_SC1_%d::tigger_restore" % target_field)


    recipe.run(["image_SC1_0_%d" % t for t in target_fields] +
               ["make_clean_mask_sc1_0_%d" % t for t in target_fields] +
               ["image_SC1_1_%d" % t for t in target_fields] +
               ["make_clean_mask_sc1_1_%d" % t for t in target_fields] +
               ["image_SC1_2_%d" % t for t in target_fields] +
               ["make_clean_mask_sc1_2_%d" % t for t in target_fields] +
               ["image_SC1_3_%d" % t for t in target_fields] +
               (["fullrestore_SC1_%d" % t for t in target_fields] if do_full_restore else []))

    ##
    # Optional : Stokes V diagnostic image - not always available so lets not
    # make this compulsory
    ##

    recipe = stimela.Recipe("Stokes V diagnostic engine", ms_dir=MSDIR)

    for target_field, target in zip(target_fields, targets):
        # Extract sources in mfs clean image to build initial sky model
        imname_prefix = cfg.obs.basename + "_STOKES_V_RESIDUE_" + plot_name[target.name]

        # it  is common belief that the Unverse is free of
        # stokes V for the most part. So any structure left
        # in it is probably calibration artifacts
        recipe.add("cab/wsclean", "wsclean_stokes_v",
                   {
                       "msname": cfg.obs.msfile,
                       "column": cfg.wsclean_v_residue.column,
                       "weight": "briggs %.2f" % (cfg.wsclean_v_residue.robust),
                       "npix": cfg.obs.im_npix,
                       "cellsize": cfg.obs.angular_resolution * cfg.obs.sampling,
                       "clean_iterations": cfg.wsclean_v_residue.clean_iterations,
                       "mgain": cfg.wsclean_v_residue.mgain,
                       "channelsout": cfg.obs.im_numchans,
                       "joinchannels": cfg.wsclean_v_residue.joinchannels,
                       "field": str(target_field),
                       "name": imname_prefix,
                       "pol": cfg.wsclean_v_residue.pol,
                   },
                   input=INPUT, output=OUTPUT,
                   label="image_stokesv_residue_%d::wsclean image STOKES "
                         "V as diagnostic" % target_field)
    try:
        recipe.run(["image_stokesv_residue_%d" % t for t in target_fields])
    except:
        vermeerkat.log.warn("Couldn't make STOKES V diagnostic image. Check "
                            "the logs for details.")
