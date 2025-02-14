from lipidRules import *


class LacCer(SphingoLipid):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in POS_ADDUCTS:
            if adduct == "[M+H]+":
                FRAGMENTS.append([PREC, 7, "precursor"])
                FRAGMENTS.append([PREC - H2O, 181, "NL H2O"])
                FRAGMENTS.append([PREC - MW("C6H12O6"), 20, "NL hexose"])

                FRAGMENTS.append([PREC - MW("C12H22O11"), 415, "NL lactose"])
                FRAGMENTS.append(
                    [PREC - MW("C12H22O11") - H2O, 99, "NL lactose + water"]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C12H22O11") - NL(self.chains[1]) + H2O,
                        151,
                        "NL lactose + acyl",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C12H22O11") - NL(self.chains[1]),
                        1000,
                        "NL lactose + acyl + H2O",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C12H22O11") - NL(self.chains[1]) + H2O - MW("CH2O"),
                        98,
                        "NL lactose + acyl + CH2O",
                    ]
                )
                FRAGMENTS.append(
                    [NL(self.chains[1]) - H2O + MW("NH3") + PROTON, 16, "fatty amide"]
                )

            if adduct == "[M+Na]+" or adduct == "[M+K]+" or adduct == ["M+NH4]+"]:
                FRAGMENTS.append([PREC, 1000, "precursor"])
                FRAGMENTS.append([PREC - MW("C6H12O6") + H2O, 46, "NL deoxyhexose"])
                FRAGMENTS.append(
                    [
                        PREC
                        - NL(self.chains[1])
                        - MW("C12H22O11")
                        - ADDUCT[adduct]
                        + PROTON,
                        18,
                        "N'': NL Lactose + acyl ketene + water + ADDUCT",
                    ]
                )
                FRAGMENTS.append(
                    [MW("C12H22O11") + ADDUCT[adduct], 12, "lactose + Adduct"]
                )
                FRAGMENTS.append(
                    [
                        MW("C12H22O11") - H2O + ADDUCT[adduct],
                        10,
                        "deoxylactose + Adduct",
                    ]
                )

        elif adduct in NEG_ADDUCTS:
            FRAGMENTS.append([PREC, 1000, "precursor"])

            if adduct == "[M+Cl]-" or adduct == "[M+FA]-" or adduct == "[M+37Cl]-":
                FRAGMENTS.append([PREC - ADDUCT[adduct] - PROTON, 1000, "pre - adduct"])
            if adduct == "[M-H]-":
                FRAGMENTS.append([PREC - MW("C6H10O5"), 145, " NL deoxyhexose"])

                FRAGMENTS.append([PREC - MW("C6H12O6"), 93, "NL hexose"])

                FRAGMENTS.append([PREC - MW("C12H20O10"), 171, "NL deoxylactose"])
                FRAGMENTS.append([PREC - MW("C12H22O11"), 1, "NL - lactose"])
                FRAGMENTS.append([MW("C12H22O11") - H2O - PROTON, 6, "deoxylactose"])

                FRAGMENTS.append([MW("C6H12O6") - PROTON, 110, "hexose"])

                FRAGMENTS.append(
                    [MW("C5H6O3") - PROTON, 167, "hexose fragment at MZ 113"]
                )
                FRAGMENTS.append(
                    [MW("C4H8O4") - PROTON, 61, "hexose fragment at MZ 119"]
                )
                FRAGMENTS.append(
                    [MW("C5H8O4") - PROTON, 38, "hexose fragment at MZ 131"]
                )
                FRAGMENTS.append(
                    [MW("C6H8O4") - PROTON, 79, "hexose fragment at MZ 143"]
                )
                FRAGMENTS.append(
                    [MW("C5H10O5") - PROTON, 13, "hexose fragment at MZ 149"]
                )
                FRAGMENTS.append(
                    [MW("C6H10O5") - PROTON, 83, "hexose fragment at MZ 161"]
                )

                FRAGMENTS.append(
                    [
                        NL(self.chains[1]) - H2O + MW("NH3") - PROTON,
                        21,
                        "sn2 fatty amide",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[1]) - H2O + MW("NH3") + MW("C2H2") - PROTON,
                        85,
                        "sn2 fatty amide + C2H2",
                    ]
                )
                FRAGMENTS.append([NL(self.chains[1]) - PROTON, 52, "sn2 fatty acid"])

            if self.chains[1][2] == 1:  # if 2-hydroxy FA
                FRAGMENTS.append(
                    [
                        NL(self.chains[1]) - H2O - MW("CO") - PROTON,
                        10,
                        "Fatty Acyl - CH2O2",
                    ]
                )  # fragment intensity not from standard
                FRAGMENTS.append(
                    [NL(self.chains[1]) - H2O + MW("H2") - PROTON, 10, "Fatty Acyl +H2"]
                )  # fragment intensity not from standard

            if self.chains[0][2] > 0:  # if not a m LCB
                FRAGMENTS.append(
                    [
                        PREC
                        - ADDUCT[adduct]
                        - NL(self.chains[1])
                        + H2O
                        - MW("C12H22O11")
                        - MW("C2H5N")
                        - PROTON,
                        14,
                        "LCB-CO2H4-CNH3",
                    ]
                )

            if self.chains[0][2] == 2:  # it a t LCB
                FRAGMENTS.append(
                    [
                        PREC
                        - ADDUCT[adduct]
                        - NL(self.chains[1])
                        + H2O
                        - MW("C12H22O11")
                        - MW("C3H7NO")
                        - PROTON,
                        1000,
                        "LCB-CO2H4-CNH3 - (CH2O if t) ",
                    ]
                )

        return FRAGMENTS


# x  = LacCer("LacCer",[[18,1,1], [17,0,0]],  adduct="[M+Na]+")
# print x.printNist()
#
# x  = LacCer("LacCer",[[18,1,1], [24,0,0]],  adduct="[M+Na]+")
# print x.printNist()
#
#
# x  = LacCer("LacCer",[[18,1,1], [24,1,0]],  adduct="[M+FA-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("LacCer", [[18,1,1], [24,1,0]], "[M+FA-H]-")
