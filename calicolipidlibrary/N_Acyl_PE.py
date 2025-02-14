from lipidRules import *


# needs work
class N_Acyl_PE(NAcylGPL):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]
        FRAGMENTS.append([PREC, 1000, "precursor"])

        if adduct in POS_ADDUCTS:
            if adduct == "[M+H]+":
                FRAGMENTS.append(
                    [PREC - MW("C2H8NO4P"), 1000, "NL ethanolamine-P"]
                )  # seems wrong
                FRAGMENTS.append(
                    [MW("C2H7NO") + ADDUCT[adduct], 1000, "ethanolamine"]
                )  # seems wrong

                for i in range(0, len(self.chains)):
                    chain = str(i + 1)
                    FRAGMENTS.append(
                        [
                            PREC - NL(self.chains[i]) + H2O,
                            1000,
                            "NL sn" + chain + "water",
                        ]
                    )
                    FRAGMENTS.append([PREC - NL(self.chains[i]), 1000, "NL sn" + chain])
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + PROTON + H2O,
                            1000,
                            "sn" + chain + " acylium ion+ water",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + PROTON,
                            1000,
                            "sn" + chain + " acylium ion",
                        ]
                    )

            if adduct == "[M+Na]+" or adduct == "[M+NH4]+" or adduct == "[M+K]+":
                FRAGMENTS.append([MW("PO4H3") + ADDUCT[adduct], 1000, "PO4H3+adduct"])
                FRAGMENTS.append(
                    [MW("C2H8NO4P") + ADDUCT[adduct], 1000, "C2H8NO4P+adduct"]
                )
                FRAGMENTS.append([PREC - MW("C2H5N"), 1000, "pre-ethanolamine"])
                # FRAGMENTS.append( [PREC-(MW("C2H8NO4P")-PROTON+ADDUCT[adduct]), 1000, "pre-141-PROTON+ADDUCT"] )
                FRAGMENTS.append([PREC - (MW("C2H8NO4P")), 1000, "pre-141"])

                for i in range(0, len(self.chains)):
                    chain = str(i + 1)
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + PROTON + H2O,
                            1000,
                            "sn" + chain + " acylium ion + water",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + PROTON,
                            1000,
                            "sn" + chain + " acylium ion",
                        ]
                    )
                    FRAGMENTS.append(
                        [PREC - NL(self.chains[i]), 1000, "NL (sn" + chain]
                    )
                    FRAGMENTS.append(
                        [
                            PREC - MW("C2H5N") - NL(self.chains[i]),
                            1000,
                            "prec-sn" + chain + "-ethanolamine-2PROTON",
                        ]
                    )  # FIX

        elif adduct in NEG_ADDUCTS:
            FRAGMENTS.append([MW("C3H9O6P") - H2O - PROTON, 1000, "glycerol3P-H2O"])
            FRAGMENTS.append([MW("H3PO4") - PROTON, 1000, "H2PO4-"])
            FRAGMENTS.append([MW("HPO3") - PROTON, 1000, "PO3-"])
            FRAGMENTS.append(
                [MW("C5H12NO5P") - PROTON, 1000, "C5H12NO5P"]
            )  # should probably remove

            if adduct == "[M+FA]-" or adduct == "[M+Cl]-":
                FRAGMENTS.append(
                    [MW("C3H9O6P") - H2O + ADDUCT[adduct], 1000, "glycerol3P-H2O"]
                )

            FRAGMENTS.append(
                [MW("C2H8NO4P") - H2O - PROTON, 1000, "Phosphoethanolamine - H20"]
            )  # should probably remove
            FRAGMENTS.append(
                [MW("C2H8NO4P") - PROTON, 1000, "Phosphoethanolamine"]
            )  # should probably remove

            FRAGMENTS.append(
                [
                    NL(self.chains[2]) - MW("O") + MW("C2H6O4PN") - PROTON,
                    1000,
                    "ethanolamide+amide",
                ]
            )
            FRAGMENTS.append(
                [
                    NL(self.chains[2])
                    - MW("O")
                    + MW("C2H6O4PN")
                    + MW("C3H4O")
                    - PROTON,
                    1000,
                    "ethanolamide+amide+glycerol",
                ]
            )
            for i in [0, 1]:
                chain = str(i + 1)
                FRAGMENTS.append(
                    [NL(self.chains[i]) - PROTON, 1000, "sn" + chain + "-RCO"]
                )
                FRAGMENTS.append(
                    [PREC - NL(self.chains[i]) + H2O, 1000, "NL sn" + chain + ""]
                )
                FRAGMENTS.append(
                    [PREC - NL(self.chains[i]), 1000, "NL sn" + chain + "-H2O"]
                )
                # should be one more which is NL of one FA and loss of amide?

        return FRAGMENTS


# x  = N_Acyl_PE("N_Acyl_PE",[[16,0,0], [18,1,0], [18,2,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("N_Acyl_PE",[[16,1,0], [18,2,0], [20,2,0]],  adduct="[M+H]+")
