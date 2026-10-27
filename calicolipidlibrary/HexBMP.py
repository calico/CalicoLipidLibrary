
from .lipidRules import *

class HexBMP(GPL):
    pos_adduct_set = ["[M+NH4]+"]
    neg_adduct_set = ["[M-H]-"]

    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in POS_ADDUCTS:
            if adduct == "[M+Na]+" or adduct == "[M+K]+" or adduct == "[M+NH4]+":

                FRAGMENTS.append([PREC, 300, "precursor"])
                FRAGMENTS.append([PREC - H2O, 50, "NL water"])
                FRAGMENTS.append([PREC - MW("C6H12O6"), 50, "precursor - hexose"])
                FRAGMENTS.append([PREC - MW("C6H10O5"), 50, "precursor - deoxyhexose"])

                FRAGMENTS.append(
                    [MW("C3H9O6P") + ADDUCT[adduct], 50, "glycerol-P + adduct"]
                )
                FRAGMENTS.append(
                    [
                        MW("C3H9O6P") - H2O + ADDUCT[adduct],
                        350,
                        "deoxyglycerol-P + Adduct",
                    ]
                )
                FRAGMENTS.append(
                    [MW("PO4H3") + ADDUCT[adduct], 50, "phosphoric acid + Adduct"]
                )

                for i in range(0, len(self.chains)):
                    chain = str(i + 1)

                    FRAGMENTS.append([NL(self.chains[i]) + MW("C3H6O2") + PROTON, 10, "Glycerol + sn" + chain])
                    FRAGMENTS.append([NL(self.chains[i]) + MW("C3H6O2") - H2O + PROTON, 1000, "Glycerol + sn" + chain + " - water"])
                    FRAGMENTS.append([NL(self.chains[i]) + MW("C3H6O2")  - H2O + ADDUCT[adduct], 10, "Glycerol + sn" + chain + " - water + Adduct"])


                    FRAGMENTS.append(
                        [
                            PREC - NL(self.chains[i]) - MW("C3H4O"),
                            10,
                            "NL (sn" + chain + " + deoxyglycerol)",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            PREC - NL(self.chains[i]) - MW("C3H6O2"),
                            10,
                            "NL (sn" + chain + " + glycerol)",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            PREC
                            - NL(self.chains[i])
                            - MW("C3H4O")
                            - MW("PO4H3")
                            - ADDUCT[adduct]
                            + PROTON,
                            200,
                            "NL (sn" + chain + " + deoxyglycerol + adduct + phosphate)",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            PREC - NL(self.chains[i]) - MW("C3H4O") - MW("PO4H3"),
                            10,
                            "NL (sn" + chain + " + deoxyglycerol + phosphate)",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            PREC - NL(self.chains[i]) - MW("C3H4O") - MW("PO4H3") + H2O,
                            10,
                            "NL (sn" + chain + " + deoxyglycerol + phosphate) + H2O",
                        ]
                    )

                    FRAGMENTS.append([PREC - NL(self.chains[i]), 2, "NL sn" + chain])
                    FRAGMENTS.append(
                        [
                            PREC - NL(self.chains[i]) + H2O,
                            10,
                            "NL sn" + chain + " ketene",
                        ]
                    )
        #These are just the exact fragments from BMPs until I have something more to work with
        elif adduct in NEG_ADDUCTS:
            FRAGMENTS.append([PREC, 501, "precursor"])

            FRAGMENTS.append([MW("C3H7O5P") - PROTON, 41, "C3H6O5P-"])

            FRAGMENTS.append([MW("C6H13O7P") - PROTON, 5, "diglycerol-phosphate"])
            FRAGMENTS.append(
                [MW("C6H11O6P") - PROTON, 1, "diglycerol-phosphate - water"]
            )

            for i in range(0, len(self.chains)):
                chain = str(i + 1)
                FRAGMENTS.append([NL(self.chains[i]) - PROTON, 500, "sn" + chain])
                FRAGMENTS.append([PREC - NL(self.chains[i]), 12, "NL (sn" + chain])
                FRAGMENTS.append(
                    [PREC - NL(self.chains[i]) + H2O, 3, "NL sn" + chain + " ketene"]
                )
                FRAGMENTS.append(
                    [
                        PREC - NL(self.chains[i]) + H2O - MW("C3H8O3"),
                        7,
                        "NL (sn" + chain + " ketene + glycerol)",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC - NL(self.chains[i]) + 2 * H2O - MW("C3H8O3"),
                        1,
                        "NL (sn" + chain + " ketene + deoxyglycerol)",
                    ]
                )

        return FRAGMENTS