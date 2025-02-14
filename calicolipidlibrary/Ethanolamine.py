from lipidRules import *


class Ethanolamine(singleAcyl):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]
        FRAGMENTS.append([PREC, 1000, "precursor"])

        if adduct in POS_ADDUCTS:  # intensites from literature, not spectra
            FRAGMENTS.append([PREC - H2O, 300, "NL water"])

            FRAGMENTS.append([PREC - MW("C2H7NO"), 300, "NL ethanolamine"])
            FRAGMENTS.append([MW("C2H7NO") + PROTON, 1000, "Ethanolamine"])

        return FRAGMENTS


# x  = Ethanolamine("Ethanolamine",[[16,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("Ethanolamine", [[16,0,0]], "[M+H]+")
