import stimela

def
    recipe = stimela.Recipe("Conversion Engine", ms_dir=MSDIR)
    # Convert
    recipe.add("cab/h5toms", "h5toms",
        {
            'hdf5files'  : [cfg.obs.h5file],
            'output-ms'  : cfg.obs.msfile,
            'model-data' : True,
            'full_pol'   : cfg.h5toms.full_pol,
        },
        input=INPUT, output=OUTPUT,
        label="convert::h5toms")

    # @mauch points out some MeerKAT h5 files contain flipped uvw
    # coordinates. Better to recalculate this from ANTENNAS
    # and be dead sure things will always work!
    recipe.add("cab/casa_fixvis", "fixvis",
        {
            "vis"       :   cfg.obs.msfile,
            "outputvis" :   cfg.obs.msfile, # into same ms please
            "reuse"     :   cfg.casa_fixvis.reuse,
        },
        input=INPUT, output=OUTPUT,
        label="recompute_uvw:: Recompute MeerKAT uvw coordinates")
    conversion_pipe = ["convert",
                       "recompute_uvw"]

    recipe.run(conversion_pipe)