from lipidRules import *


class Alkyl_LPS(AlkylLysoGPL):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in POS_ADDUCTS:
            FRAGMENTS.append([PREC, 1000, "precursor"])

        elif adduct in NEG_ADDUCTS:
            FRAGMENTS.append([PREC, 1000, "precursor"])

            FRAGMENTS.append([MW("H3PO4") - PROTON, 1000, "Phosphate"])
            FRAGMENTS.append([MW("HPO3") - PROTON, 1000, "phosphite"])
            FRAGMENTS.append([PREC - MW("C3H5NO2"), 1000, "NL Serine"])

            FRAGMENTS.append(
                [MW("C3H7O5P") - PROTON, 1000, "glycerol phosphate - water"]
            )
            FRAGMENTS.append([MW("C3H9O6P") - PROTON, 1000, "glycerol phosphate"])

        return FRAGMENTS


# x  = Alkyl_LPS("Alkyl_LPS",[[16,0,0]],  adduct="[M+H]+")
# print x.printNist()

# calicolipidlibrary.print_spectrum("Alkyl_LPS",[[16,0,0]], "[M+H]+")
