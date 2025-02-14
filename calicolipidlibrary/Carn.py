from lipidRules import *


class Carn(singleAcyl):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in POS_ADDUCTS:
            FRAGMENTS.append([PREC, 837, "precursor"])
            FRAGMENTS.append([MW("C3H9N") + PROTON, 83, "trimethylamine"])

            FRAGMENTS.append([PREC - MW("C3H9N"), 129, "NL trimethylamine"])

            FRAGMENTS.append([MW("C4H4O2") + ADDUCT[adduct], 890, "C4H5O2 + adduct "])
            FRAGMENTS.append([MW("C7H13NO2") + ADDUCT[adduct], 32, "deoxycarnitine"])
            FRAGMENTS.append(
                [NL(self.chains[0]) - H2O + ADDUCT[adduct], 40, "acylium ion "]
            )

        return FRAGMENTS


# x  = Carn("Carn",[[12,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
#
#
# x  = Carn("Carn",[[16,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("Carn", [[16,0,0]], "[M+H]+")
