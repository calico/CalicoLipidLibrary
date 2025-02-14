from lipidRules import *


class PS(GPL):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in POS_ADDUCTS:
            if adduct == "[M+H]+":
                FRAGMENTS.append([PREC, 1, "precursor"])

                FRAGMENTS.append([PREC - H2O, 3, "NL water"])  # C3H8NO6P
                FRAGMENTS.append(
                    [PREC - MW("C3H8NO6P"), 1000, "NL phosphoserine"]
                )  # C3H8NO6P
                FRAGMENTS.append([MW("C3H7NO3") + ADDUCT[adduct], 3, "Serine"])  # Ser
                FRAGMENTS.append(
                    [MW("C3H7NO3") + ADDUCT[adduct] - H2O, 8, "Serine - water"]
                )  # Ser
                FRAGMENTS.append(
                    [MW("C3H9O6P") + ADDUCT[adduct] - H2O, 6, "Glycerol-P - water"]
                )

                for i in range(0, len(self.chains)):
                    chain = str(i + 1)
                    FRAGMENTS.append(
                        [
                            PREC - NL(self.chains[i]) + H2O,
                            2,
                            "NL sn" + chain + "acylium",
                        ]
                    )
                    FRAGMENTS.append([PREC - NL(self.chains[i]), 2, "NL sn" + chain])
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O + PROTON,
                            26,
                            "sn" + chain + " acylium ion",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O - H2O + PROTON,
                            2,
                            "sn" + chain + " acylium ion - H2O",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O + MW("C3H6O2") + ADDUCT[adduct],
                            24,
                            "sn" + chain + " + glycerol",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O + MW("C3H5NO2") + ADDUCT[adduct],
                            23,
                            "sn" + chain + " + serine",
                        ]
                    )

            if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
                FRAGMENTS.append([PREC, 357, "precursor"])
                FRAGMENTS.append(
                    [
                        MW("C3H7NO3") + ADDUCT[adduct] - H2O,
                        43,
                        "Serine - water + adduct",
                    ]
                )  # Ser
                FRAGMENTS.append(
                    [
                        PREC - MW("C3H8NO6P") + PROTON - ADDUCT[adduct],
                        202,
                        "NL P-Ser + adduct",
                    ]
                )
                FRAGMENTS.append([PREC - MW("C3H8NO6P"), 59, "NL phosphoSer"])
                FRAGMENTS.append(
                    [MW("C3H8NO6P") + ADDUCT[adduct], 1000, "P-Ser+adduct"]
                )
                FRAGMENTS.append([MW("H3PO4") + ADDUCT[adduct], 164, " H3PO4+adduct"])

                for i in range(0, len(self.chains)):
                    chain = str(i + 1)
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O + PROTON,
                            8,
                            "sn" + chain + " acylium ion",
                        ]
                    )  # keep
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O - H2O + PROTON,
                            1,
                            "sn" + chain + " acylium ion - H2O",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            PREC - NL(self.chains[i]) - MW("C3H7NO3") + H2O,
                            7,
                            "NL sn" + chain + " deoxyserine",
                        ]
                    )  # keep

        elif adduct in NEG_ADDUCTS:
            FRAGMENTS.append([PREC, 261, "precursor"])

            FRAGMENTS.append([MW("H3PO4") - PROTON, 22, "H2PO4-"])
            FRAGMENTS.append([PREC - MW("C3H5NO2"), 477, "NL Serine"])
            FRAGMENTS.append([MW("C3H7O5P") - H2O - PROTON, 5, "glycerol-P - 2xwater"])
            FRAGMENTS.append([MW("C3H7O5P") - PROTON, 290, "glycerol-P - water"])
            FRAGMENTS.append([MW("C3H9O6P") - PROTON, 5, "Glycerol-P"])

            for i in range(0, len(self.chains)):
                chain = str(i + 1)
                FRAGMENTS.append(
                    [NL(self.chains[i]) - PROTON, 500, "sn" + chain + " RCO2-"]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C3H5NO2") - NL(self.chains[i]) + H2O,
                        58,
                        "NL (Serine + sn" + chain + ")",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C3H5NO2") - NL(self.chains[i]),
                        160,
                        "NL Serine + sn" + chain + " + H2O)",
                    ]
                )

        return FRAGMENTS


# x  = PS("PS",[[16,0,0], [18,1,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("PS",[[16,0,0], [18,1,0]],  adduct="[M+H]+")
