from lipidRules import *


class DG(GPL):
    pos_adduct_set = ["[M+Na]+","[M+NH4]+"]

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
                            NL(self.chains[i]) + PROTON - H2O + MW("C3H6O2"),
                            14,
                            "sn" + chain + " + deoxyglycerol",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + ADDUCT[adduct] - H2O + MW("C3H6O2"),
                            5,
                            "sn" + chain + " + deoxyglycerol + adduct",
                        ]
                    )

        return FRAGMENTS


# x  = DG("DG",[[16,0,0], [20,1,0]],  adduct="[M+NH4]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("DG", [[16,0,0], [20,1,0]], "[M+NH4]+")
