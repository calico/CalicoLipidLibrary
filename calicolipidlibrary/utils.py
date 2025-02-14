from AcGM2 import *
from AcGM3 import *
from Alkyl_LPC import *
from Alkyl_LPE import *
from Alkyl_LPS import *
from Alkyl_PC import *
from Alkyl_PE import *
from Alkyl_PS import *
from GB3 import *
from BDP import *
from BMP import *
from Carn import *
from CDP_DG import *
from CE import *
from Ceramide import *
from Ceramide_P import *
from CL import *
from CPE import *
from CPI import *
from DG import *
from DGDG import *
from DGTS import *
from DMPE import *
from ErgE import *
from Ethanolamine import *
from FA import *
from FAHFA import *
from GcGM2 import *
from GcGM3 import *
from HemiBMP import *
from HexCer import *
from LacCer import *
from LCB_P import *
from LCB import *
from LPA import *
from LPC import *
from LPE import *
from LPG import *
from LPI import *
from LPS import *
from LysoCL import *
from LysoCPE import *
from LysoCPI import *
from LysoHexCer import *
from LysoSM import *
from MG import *
from MGDG import *
from MIP2C import *
from MIPC import *
from MMPE import *
from N_Acyl_PE import *
from N_Acyl_PS import *
from PA import *
from PC import *
from PE import *
from PG import *
from PI import *
from PS import *
from SM import *
from Sulfatide import *
from Taurine import *
from TG import *

from lipidRules import *


def print_spectrum(lipid_class, chains=[], adduct=""):
    if adduct not in ALL_ADDUCTS:
        print("Adduct '" + adduct + "' is not a supported adduct form.")
        print("Valid adducts are:")
        for adduct in ALL_ADDUCTS:
            print(adduct)
        return

    if lipid_class not in ALL_LIPID_CLASSES:
        print("Lipid Class '" + lipid_class + "' is not a supported lipid class.")
        print("Suported lipid classes are:")
        for lipid in ALL_LIPID_CLASSES:
            print(lipid)
        return

    eval_str = (
        lipid_class
        + "('"
        + lipid_class
        + "', "
        + chains.__str__()
        + ", '"
        + adduct
        + "').printNist()"
    )

    print(eval(eval_str))
