from lipidRules import *


class LPA(LysoGPL):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in NEG_ADDUCTS:
            if self.chains[0][0] != 0:
                i = 0
            else:
                i = 1
            chain = str(i + 1)

            FRAGMENTS.append([PREC, 54, "precursor"])

            FRAGMENTS.append([MW("H3PO4") - PROTON, 20, "phosphate"])
            FRAGMENTS.append([MW("HPO3") - PROTON, 223, "phosphite"])
            FRAGMENTS.append(
                [MW("C3H5O5P") - PROTON, 6, "glycerol phosphate - water - H2"]
            )
            FRAGMENTS.append(
                [MW("C3H7O5P") - PROTON, 1000, "glycerol phosphate - water"]
            )
            FRAGMENTS.append([MW("C3H9O6P") - PROTON, 32, "glycerol phosphate"])

            FRAGMENTS.append([NL(self.chains[i]) - PROTON, 17, "sn" + chain])

        return FRAGMENTS


# x  = LPA("LPA",[[16,0,0],[0,0,0]],  adduct="[M-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("LPA", [[16,0,0],[0,0,0]], "[M-H]-")
