import os
from scipy.optimize import curve_fit

import numpy as np

__CALIBRATOR_DB = None

def calibrator_database():
    """ Return the Southern standard calibrator database """

    global __CALIBRATOR_DB

    # Do a lazy load
    if __CALIBRATOR_DB is not None:
        return __CALIBRATOR_DB

    # OK its not loaded, read it in
    import vermeerkat
    import vermeerkat.dispatch_crew.caltable_parser as vmcp

    # There isn't a Southern standard in CASA
    # so construct a little database of them for reference
    vermeerkat.log.info("Parsing calibrator table")

    ref_table = os.path.join(vermeerkat.install_path(), "data/southern_calibrators.txt")
    __CALIBRATOR_DB = vmcp.read_caltable(ref_table)

    vermeerkat.log.info("Found the following reference calibrators (in GHz format):")
    vermeerkat.log.info(vmcp.format_calibrator_db(__CALIBRATOR_DB))

    return __CALIBRATOR_DB

def convert_pb_to_casaspi(vlower, vupper, v0, a, b, c, d):
    """
    Coverts between the different conventions:
    PB: 10 ** [a + b * log10(v) + c * log10(v) ** 2 + d * log10(v) ** 3]
    CASA/Meqtrees SPI: S(v0) * (v/v0) ** [a' + b'*log10(v/v0) + c'*log10(v/v0) ** 2 + d'*log10(v/v0) ** 3]

    args:
    :vlower, vupper: range (same unit as a, b, c, d coefficients!) to fit for a',b',c',d'
    :v0: reference frequency (same unit as vlower, vupper!)
    :a,b,c,d: PB coefficients (for the unit used in vlower, vupper and v0!)
    """
    if vlower > vupper:
        raise ValueError("vlower must be lower than vupper")

    def pbspi(v, a, b, c, d):
        return 10 ** (a + b * np.log10(v) + c * np.log10(v) ** 2 + d * np.log10(v) ** 3)
    def casaspi(v, v0, I, a, b, c, d):
        return I * (v/v0) ** (a + b * np.log10(v/v0) + c * np.log10(v/v0) ** 2 + d * np.log10(v/v0) ** 3)

    I = pbspi(v0, a, b, c, d)

    v = np.linspace(vlower, vupper, 10000)
    popt, pcov = curve_fit(lambda v, a, b, c, d: casaspi(v, v0, I, a, b, c, d), v, pbspi(v,a,b,c,d))
    perr = np.sqrt(np.diag(pcov))
    assert np.all(perr < 1.0e-6)

    # returns (S(v0), a', b', c', d')
    return I, popt[0], popt[1], popt[2], popt[3]
