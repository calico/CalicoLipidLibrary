from lipidRules import *


class DMPE(GPL):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in POS_ADDUCTS:
            if adduct == "[M+H]+":
                FRAGMENTS.append([PREC, 1, "precursor"])
                FRAGMENTS.append(
                    [MW("C4H12NO4P") + ADDUCT[adduct], 443, "DMPE phosphate"]
                )
                FRAGMENTS.append([MW("C4H9N") + PROTON, 495, "DMPE - water"])
                FRAGMENTS.append([PREC - MW("C4H12NO4P"), 1000, "pre-MMPE-Phosphate"])

                i = 0
                chain = str(i + 1)
                # FRAGMENTS.append( [PREC - NL(self.chains[i])-H2O, 1000, "NL sn" + chain + " + water"])
                # FRAGMENTS.append( [PREC - NL(self.chains[i]), 1000, "NL sn" + chain ])
                FRAGMENTS.append(
                    [
                        NL(self.chains[i]) + PROTON - H2O,
                        27,
                        "sn" + chain + " acylium ion ",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[i]) + PROTON - H2O - H2O,
                        4,
                        "sn" + chain + " acylium ion - water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[i]) + PROTON + MW("C3H4O"),
                        5,
                        "sn" + chain + " deoxyglycerol",
                    ]
                )

                i = 1
                chain = str(i + 1)
                # FRAGMENTS.append( [PREC - NL(self.chains[i])-H2O, 1000, "NL sn" + chain + " + water"])
                # FRAGMENTS.append( [PREC - NL(self.chains[i]), 1000, "NL sn" + chain ])
                FRAGMENTS.append(
                    [
                        NL(self.chains[i]) + PROTON - H2O,
                        38,
                        "sn" + chain + " acylium ion ",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[i]) + PROTON - H2O - H2O,
                        15,
                        "sn" + chain + " acylium ion - water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[i]) + PROTON + MW("C3H4O"),
                        5,
                        "sn" + chain + " deoxyglycerol",
                    ]
                )

            if adduct == "[M+Na]+" or adduct == "[M+NH4]+" or adduct == "[M+K]+":
                FRAGMENTS.append([PREC, 1000, "precursor"])
                FRAGMENTS.append(
                    [MW("C4H12NO4P") + ADDUCT[adduct], 1000, "DMPE phosphate"]
                )
                FRAGMENTS.append([MW("C4H9N") + PROTON, 1000, "DMPE - water"])
                FRAGMENTS.append(
                    [MW("PO4H3") + ADDUCT[adduct], 1000, "phosphoric acid + adduct"]
                )
                FRAGMENTS.append([PREC - MW("C4H9N"), 1000, "NL DMPE"])
                FRAGMENTS.append(
                    [PREC - (MW("C4H12NO4P")), 1000, "NL DMPE + phosphate"]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C4H12NO4P") - ADDUCT[self.adduct] + PROTON,
                        1000,
                        "NL DMPE + phosphate + ADDUCT",
                    ]
                )

                for i in range(0, len(self.chains)):
                    chain = str(i + 1)
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + PROTON - H2O - H2O,
                            1000,
                            "sn" + chain + " acylium ion - water",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + PROTON - H2O,
                            1000,
                            "sn" + chain + " acylium ion",
                        ]
                    )
                    FRAGMENTS.append([PREC - NL(self.chains[i]), 1000, "NL sn" + chain])
                    FRAGMENTS.append(
                        [
                            PREC - MW("C4H9N") - NL(self.chains[i]),
                            1000,
                            "NL sn" + chain + " + DMPE",
                        ]
                    )

        elif adduct in NEG_ADDUCTS:
            FRAGMENTS.append([PREC, 1000, "precursor"])
            # these look just lie PC's that demethylate in source

        return FRAGMENTS


# x  = DMPE("DMPE",[[16,0,0], [16,1,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("DMPE", [[16,0,0], [16,1,0]], "[M+H]+")
