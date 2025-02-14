from lipidRules import *


class ErgE(singleAcyl):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in POS_ADDUCTS:
            FRAGMENTS.append([MW("C27H44") + PROTON, 1000, "Ergosterol - water"])

            if adduct == "[M+Na]+" or adduct == "[M+K]+" or adduct == "[M+NH4]+":
                FRAGMENTS.append([PREC, 1000, "precursor"])
                FRAGMENTS.append(
                    [MW("C28H44O") - H2O + PROTON, 8, "Ergosterol - water"]
                )
                FRAGMENTS.append(
                    [NL(self.chains[0]) + ADDUCT[adduct], 249, "Fatty acyl + adduct"]
                )
                # FRAGMENTS.append([MW("C11H14")+PROTON, 8, "C11H15+"])
                if self.chains[0][2] > 0:
                    FRAGMENTS.append(
                        [
                            NL(self.chains[0]) + ADDUCT[adduct] - H2O,
                            87,
                            "Fatty acyl + adduct - H2O",
                        ]
                    )

        return FRAGMENTS


# x  = ErgE("ErgE",[[16,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("ErgE", [[16,0,0]], "[M+H]+")
