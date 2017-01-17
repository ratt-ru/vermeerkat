import re
import numpy as np

def read_caltable(filename):
    """
        Read calibrator database (specified in MHz)
        and returns a dictionary containing the following
        :filename: filename of caltable database
        :returns: for every source (name = key):
                    Epoch, RA, Declination, 
                    a_ghz, b_ghz, c_ghz, d_ghz,
                    a_mhz, b_mhz, c_mhz, d_mhz
    """
    calibrator_db = {}
    with open(filename) as f:
        line=f.readline()
        ln_no = 1
        while line:
            #discard comments
            command = line.split("//")[0]

            #empty line ?
            if command.strip() == "":
                line = f.readline()
                ln_no += 1
                continue

            #source ?
            valset = re.match(r"^name=(?P<name>[0-9A-Za-z\-+_]+)[ ]+"
                              r"epoch=(?P<epoch>[0-9]+)[ ]+"
                              r"ra=(?P<ra>[+\-]?[0-9]+h[0-9]+m[0-9]+(?:.[0-9]+)?s)[ ]+"
                              r"dec=(?P<decl>[+\-]?[0-9]+d[0-9]+m[0-9]+(?:.[0-9]+)?s)[ ]+"
                              r"a=(?P<a>[+\-]?[0-9]+(?:.[0-9]+)?)[ ]+"
                              r"b=(?P<b>[+\-]?[0-9]+(?:.[0-9]+)?)[ ]+"
                              r"c=(?P<c>[+\-]?[0-9]+(?:.[0-9]+)?)[ ]+"
                              r"d=(?P<d>[+\-]?[0-9]+(?:.[0-9]+)?)$", 
                              command)
            #else illegal
            if not valset:
                raise RuntimeError("Illegal line encountered while parsing"
                                   "southern standard at line %d:'%s'" %
                                   (ln_no, line))

            # parse sources (spectra in MHz) 
            name = valset.group("name")
            epoch = int(valset.group("epoch"))
            ra = valset.group("ra")
            valset_ra = re.match(r"^(?P<h>[+\-]?[0-9]+)h"
                                 r"(?P<m>[0-9]+)m"
                                 r"(?P<s>[0-9]+(?:.[0-9]+)?)s$",
                                 ra)
            ra = np.deg2rad((float(valset_ra.group("h")) +
                             float(valset_ra.group("m")) / 60.0 +
                             float(valset_ra.group("s")) / 3600) / 24.0 * 360)
            decl = valset.group("decl")
            valset_decl = re.match(r"^(?P<d>[+\-]?[0-9]+)d"
                                   r"(?P<m>[0-9]+)m"
                                   r"(?P<s>[0-9]+(?:.[0-9]+)?)s$",
                                   decl)
            decl = np.deg2rad(float(valset_decl.group("d")) + \
                              float(valset_decl.group("m")) + \
                              float(valset_decl.group("s")))

            a = float(valset.group("a"))
            b = float(valset.group("b"))
            c = float(valset.group("c"))
            d = float(valset.group("d"))

            # convert models to Perley Butler GHz format
            k = np.log10(1000)
            ag = a + (b * k) + (c * k ** 2) + (d * k ** 3)
            bg = b + (2 * c * k) + (3 * d * k ** 2)
            cg = c + (3 * d * k)
            dg = d

            calibrator_db[name] = {"epoch": epoch,
                                   "ra": ra,
                                   "decl": decl,
                                   "a_ghz": ag,
                                   "b_ghz": bg,
                                   "c_ghz": cg,
                                   "d_ghz": dg,
                                   "a_mhz": a,
                                   "b_mhz": b,
                                   "c_mhz": c,
                                   "d_mhz": d}
            # finally parse next line
            line = f.readline()
            ln_no += 1

        return calibrator_db
