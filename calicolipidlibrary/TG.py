import sys
from lipidRules import *


class TG(Triglyceride):
    pos_adduct_set = ["[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+Li]+"]

    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in POS_ADDUCTS:
            if adduct == "[M+H]+":
                FRAGMENTS.append([PREC, 1000, "precursor"])

            if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
                FRAGMENTS.append([PREC, 1000, "precursor"])

                for i in range(0, len(self.chains)):
                    chain = str(i + 1)
                    FRAGMENTS.append(
                        [
                            PREC - NL(self.chains[i]) - ADDUCT[adduct] + PROTON,
                            61,
                            "NL sn" + chain + " + adduct",
                        ]
                    )
                    FRAGMENTS.append([PREC - NL(self.chains[i]), 56, "NL sn" + chain])
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + PROTON - H2O,
                            10,
                            "sn" + chain + " acylium ion",
                        ]
                    )

        return FRAGMENTS


# x  = TG("TG",[[15,0,0], [18,1,0], [15,0,0]],  adduct="[M+Na]+")
# print x.printNist()
#
# x  = TG("TG",[[15,0,0], [18,1,0], [15,0,0]],  adduct="[M+NH4]+")
# print x.printNist()
#
# x  = TG("TG",[[15,0,0], [18,1,0], [15,0,0]],  adduct="[M+Li]+")
# print x.printNist()

# calicolipidlibrary.print_spectrum("TG",[[16,0,0], [18,1,0], [18,0,0]], "[M+NH4]+")
