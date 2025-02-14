from lipidRules import *


class HemiBMP(LysoCardioLipin):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in NEG_ADDUCTS:
            FRAGMENTS.append([PREC, 447, "precursor"])

            FRAGMENTS.append([MW("C3H7O5P") - PROTON, 41, "C3H6O5P-"])

            for i in range(0, len(self.chains)):
                chain = str(i + 1)
                FRAGMENTS.append(
                    [NL(self.chains[i]) - PROTON, 333, "sn" + chain + " RCO2-"]
                )
                FRAGMENTS.append(
                    [PREC - NL(self.chains[i]), 3, "NL (sn" + chain + " + water)"]
                )
                FRAGMENTS.append([PREC - NL(self.chains[i]) + H2O, 4, "NL sn" + chain])

                FRAGMENTS.append(
                    [
                        NL(self.chains[i]) + MW("C3H7O5P") - PROTON,
                        1,
                        "sn" + chain + " LPA",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[i]) + MW("C3H7O5P") - H2O - PROTON,
                        8,
                        "sn" + chain + " LPA - water",
                    ]
                )

        return FRAGMENTS


# x  = HemiBMP("HemiBMP",[[18,0,0], [22,6,0], [2,0,0]],  adduct="[M-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("HemiBMP", [[18,0,0], [22,6,0], [2,0,0]], "[M-H]-")
