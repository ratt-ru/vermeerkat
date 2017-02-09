import os

__CALIBRATOR_DB = None

def calibrator_database():
    """ Return the Southern standard calibrator database """

    global __CALIBRATOR_DB

    # Do a lazy load
    if __CALIBRATOR_DB is not None:
        return __CALIBRATOR_DB

    # OK its not loaded, read it in
    import vermeerkat
    import vermeerkat.caltable_parser as vmcp

    # There isn't a Southern standard in CASA
    # so construct a little database of them for reference
    vermeerkat.log.info("Parsing calibrator table")

    ref_table = os.path.join(vermeerkat.install_path(), "southern_calibrators.txt")
    __CALIBRATOR_DB = vmcp.read_caltable(ref_table)

    vermeerkat.log.info("Found the following reference calibrators (in GHz format):")
    vermeerkat.log.info(vmcp.format_calibrator_db(__CALIBRATOR_DB))

    return __CALIBRATOR_DB