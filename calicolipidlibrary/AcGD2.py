from lipidRules import *


class AcGD2(SphingoLipid):
    neg_adduct_set = ["[M-H]-", "[M-2H]2-"]
    pos_adduct_set = ["[M+H]+", "[M+Na]+"]

    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = (MASS + ADDUCT[self.adduct]) / abs(ADDUCT_CHARGE[self.adduct])

        if adduct in POS_ADDUCTS:
            if adduct == "[M+H]+":
                FRAGMENTS.append([PREC, 1, "precursor"])
                FRAGMENTS.append([PREC - H2O, 15, "NL water"])
                FRAGMENTS.append([PREC - MW("C8H15NO6"), 12, "NL N-Ac-Galactosamine"])
                FRAGMENTS.append(
                    [PREC - MW("C11H19NO9"), 21, "NL N-Ac-Neuraminic Acid"]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C11H19NO9") - MW("C8H15NO6") + H2O,
                        23,
                        "NL deoxy-N-Ac-Neuraminic Acid + N-Ac-Galactosamine ",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC - 2 * MW("C11H19NO9") + H2O,
                        119,
                        "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid",
                    ]
                )

                # FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+2*H2O, 98, "NL deoxy-N-Ac-Neuraminic Acid + deoxyGalactose + N-Ac-Galactosamine "] )
                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - MW("C6H12O6")
                        + 3 * H2O,
                        78,
                        "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H12O6")
                        + 3 * H2O,
                        51,
                        "Ceramide - 2x water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H12O6")
                        + 4 * H2O,
                        230,
                        "Ceramide - water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H12O6")
                        + 5 * H2O,
                        2,
                        "Ceramide",
                    ]
                )

                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H12O6")
                        + 4 * H2O
                        - NL(self.chains[1])
                        + H2O,
                        59,
                        "sphingoid base ",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H12O6")
                        + 3 * H2O
                        - NL(self.chains[1])
                        + H2O,
                        542,
                        "sphingoid base - water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H12O6")
                        + 4 * H2O
                        - NL(self.chains[1])
                        + H2O
                        - MW("CH2O"),
                        45,
                        "sphingoid base - CH2O",
                    ]
                )

                FRAGMENTS.append(
                    [
                        NL(self.chains[1]) - H2O + MW("NH3") + PROTON,
                        11,
                        "sn2 fatty amide",
                    ]
                )

                FRAGMENTS.append([MW("C11H19NO9") + PROTON, 18, "N-Ac-Neuraminic Acid"])
                FRAGMENTS.append(
                    [MW("C11H19NO9") - H2O + PROTON, 583, "deoxy-N-Ac-Neuraminic Acid"]
                )
                FRAGMENTS.append(
                    [
                        MW("C11H19NO9") - 2 * H2O + PROTON,
                        1000,
                        "deoxy-N-Ac-Neuraminic Acid - water ",
                    ]
                )
                FRAGMENTS.append(
                    [MW("C8H15NO6") - H2O + PROTON, 488, "deoxy-N-Ac-glucosamine"]
                )
                FRAGMENTS.append(
                    [
                        MW("C8H15NO6") - 2 * H2O + PROTON,
                        204,
                        "deoxy-N-Ac-glucosamine - water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        MW("C6H12O6") + MW("C11H19NO9") - 2 * H2O + PROTON,
                        27,
                        "deoxygalactose + deoxy-N-Ac-Neuraminic Acid",
                    ]
                )
                # FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3*H2O + PROTON, 16, "deoxy-N-Ac-galcatosamine + deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
                FRAGMENTS.append(
                    [
                        MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + PROTON,
                        73,
                        "deoxy-N-Ac-glucosamine + deoxy-galactose",
                    ]
                )
                FRAGMENTS.append(
                    [
                        MW("C8H15NO6") + 2 * MW("C6H12O6") - 3 * H2O + PROTON,
                        2,
                        "deoxy-N-Ac-glucosamine + deoxy-galactose + deoxy-glucose",
                    ]
                )

        # omitted [M+2H]2+ as the reference specctrum was poor quality

        if adduct == "[M+Na]+" or adduct == "[M+K]+":
            FRAGMENTS.append([PREC, 66, "precursor"])
            FRAGMENTS.append([PREC - MW("CO2"), 3, "NL CO2"])
            FRAGMENTS.append([PREC - H2O, 3, "NL water"])
            FRAGMENTS.append([PREC - MW("C8H15NO6") + H2O, 2, "NL N-Ac-Galactosamine"])

            FRAGMENTS.append(
                [PREC - MW("C11H19NO9") + H2O, 66, "NL deoxy-N-Ac-Neuraminic Acid"]
            )
            FRAGMENTS.append([PREC - MW("C11H19NO9"), 2, "NL N-Ac-Neuraminic Acid"])
            FRAGMENTS.append(
                [
                    PREC - 2 * MW("C11H19NO9") + 2 * H2O,
                    1000,
                    "NL 2x deoxy-N-Ac-Neuraminic Acid",
                ]
            )
            FRAGMENTS.append(
                [
                    PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") + 3 * H2O,
                    257,
                    "NL 2x deoxy-N-Ac-Neuraminic Acid + deoxy-N-Ac-Galactosamine",
                ]
            )
            FRAGMENTS.append(
                [
                    PREC
                    - 2 * MW("C11H19NO9")
                    - MW("C8H15NO6")
                    - MW("C6H12O6")
                    + 4 * H2O,
                    40,
                    "NL 2x deoxy-N-Ac-Neuraminic Acid + deoxy-Galactose + deoxy-N-Ac-Galactosamine",
                ]
            )
            FRAGMENTS.append(
                [
                    PREC
                    - 2 * MW("C11H19NO9")
                    - MW("C8H15NO6")
                    - 2 * MW("C6H12O6")
                    + 5 * H2O,
                    2,
                    "NL 2x deoxy-N-Ac-Neuraminic Acid + 2x deoxyGalactose + deoxy-N-Ac-Galactosamine",
                ]
            )

            FRAGMENTS.append(
                [
                    MW("C8H15NO6") - H2O + ADDUCT[adduct],
                    6,
                    "deoxy-N-Ac-galactosamine + Adduct",
                ]
            )
            FRAGMENTS.append(
                [MW("C11H19NO9") - H2O + PROTON, 1, "deoxy-N-Ac-Neuraminic Acid"]
            )
            FRAGMENTS.append(
                [
                    MW("C11H19NO9") - H2O + ADDUCT[adduct],
                    12,
                    "deoxy-N-Ac-Neuraminic Acid + Adduct",
                ]
            )
            FRAGMENTS.append(
                [
                    MW("C11H19NO9") - 2 * H2O + PROTON,
                    2,
                    "deoxy-N-Ac-Neuraminic Acid - water ",
                ]
            )

            FRAGMENTS.append(
                [
                    MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + PROTON,
                    12,
                    "deoxy-N-Ac-galactosamine + deoxygalactose",
                ]
            )
            FRAGMENTS.append(
                [
                    MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + ADDUCT[adduct],
                    4,
                    "deoxy-N-Ac-galactosamine + deoxygalactose + ADDUCT",
                ]
            )

        elif adduct in NEG_ADDUCTS:
            if adduct == "[M+Cl]-" or adduct == "[M+FA]-":
                FRAGMENTS.append([PREC - ADDUCT[adduct] - PROTON, 1000, "NL adduct"])
            if adduct == "[M-H]-":
                FRAGMENTS.append([PREC, 179, "precursor"])
                FRAGMENTS.append([PREC - H2O, 3, "NL water"])
                FRAGMENTS.append([PREC - MW("CO2"), 25, "NL CO2"])
                FRAGMENTS.append(
                    [
                        NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON,
                        8,
                        "N-acyl-ethanolamine - H2O",
                    ]
                )

                FRAGMENTS.append([MW("C11H19NO9") - PROTON, 6, "N-Ac-Neuraminic Acid"])
                FRAGMENTS.append(
                    [MW("C11H19NO9") - H2O - PROTON, 481, "deoxy-N-Ac-Neuraminic Acid"]
                )
                FRAGMENTS.append(
                    [
                        MW("C11H19NO9") + MW("C11H19NO9") - 2 * H2O - PROTON,
                        146,
                        "2x deoxy-N-Ac-Neuraminic Acid",
                    ]
                )
                FRAGMENTS.append(
                    [
                        MW("C11H19NO9")
                        + MW("C11H19NO9")
                        - MW("CO2")
                        - 2 * H2O
                        - PROTON,
                        56,
                        "2x deoxy-N-Ac-Neuraminic Acid - CO2",
                    ]
                )
                FRAGMENTS.append(
                    [
                        MW("C11H19NO9")
                        + MW("C11H19NO9")
                        - MW("CO2")
                        - 3 * H2O
                        - PROTON,
                        6,
                        "2x deoxy-N-Ac-Neuraminic Acid - CO2 - water",
                    ]
                )

                FRAGMENTS.append(
                    [
                        PREC - MW("C11H19NO9") + H2O,
                        1000,
                        "NL deoxy-N-Ac-Neuraminic Acid",
                    ]
                )
                # FRAGMENTS.append( [PREC-MW("C11H19NO9"), 4, "NL N-Ac-Neuraminic Acid"] )
                # FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("CO2")+H2O, 2, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )
                FRAGMENTS.append(
                    [
                        PREC - 2 * MW("C11H19NO9") + 2 * H2O,
                        162,
                        "NL 2x deoxy-N-Ac-Neuraminic Acid",
                    ]
                )

                # FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 2, "Fatty Amide - NH3 - H"])
                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        + 3 * H2O,
                        21,
                        "Ceramide",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - MW("C6H10O5")
                        + 2 * H2O,
                        2,
                        "Glucosyl ceramide - water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - MW("C6H10O5")
                        + 3 * H2O,
                        14,
                        "Glucosyl ceramide",
                    ]
                )
                # FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") -  MW("C6H10O5") + 3 * H2O, 48, "Lactosyl ceramide"])

            if adduct == "[M-2H]2-":
                SINGLE = MASS - PROTON
                FRAGMENTS.append([PREC, 1000, "precursor"])
                FRAGMENTS.append([PREC - H2O / 2, 2, "NL water 2-"])
                FRAGMENTS.append([PREC - MW("CO2") / 2, 14, "NL CO2 2-"])

                FRAGMENTS.append(
                    [
                        NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON,
                        16,
                        "N-acyl-ethanolamine - H2O",
                    ]
                )
                FRAGMENTS.append([MW("C11H19NO9") - PROTON, 16, "N-Ac-Neuraminic Acid"])
                FRAGMENTS.append(
                    [MW("C11H19NO9") - H2O - PROTON, 488, "deoxy-N-Ac-Neuraminic Acid"]
                )
                FRAGMENTS.append(
                    [
                        MW("C11H19NO9") + MW("C11H19NO9") - 2 * H2O - PROTON,
                        223,
                        "2x deoxy-N-Ac-Neuraminic Acid",
                    ]
                )
                FRAGMENTS.append(
                    [
                        MW("C11H19NO9")
                        + MW("C11H19NO9")
                        - MW("CO2")
                        - 2 * H2O
                        - PROTON,
                        11,
                        "2x deoxy-N-Ac-Neuraminic Acid - CO2",
                    ]
                )
                FRAGMENTS.append(
                    [
                        MW("C11H19NO9")
                        + MW("C11H19NO9")
                        - MW("CO2")
                        - 3 * H2O
                        - PROTON,
                        4,
                        "2x deoxy-N-Ac-Neuraminic Acid - CO2 - water",
                    ]
                )

                FRAGMENTS.append(
                    [
                        SINGLE - MW("C11H19NO9") + H2O,
                        99,
                        "NL deoxy-N-Ac-Neuraminic Acid",
                    ]
                )

                # FRAGMENTS.append( [SINGLE-MW("C11H19NO9"), 1, "NL N-Ac-Neuraminic Acid"] )
                # FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("CO2")+H2O, 2, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )
                FRAGMENTS.append(
                    [
                        SINGLE - 2 * MW("C11H19NO9") + 2 * H2O,
                        222,
                        "NL 2x deoxy-N-Ac-Neuraminic Acid",
                    ]
                )

                # FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 2, "Fatty Amide - NH3 - H"])
                FRAGMENTS.append(
                    [
                        SINGLE
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        + 3 * H2O,
                        41,
                        "Ceramide",
                    ]
                )
                FRAGMENTS.append(
                    [
                        SINGLE
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        + 2 * H2O,
                        2,
                        "Ceramide - H2O",
                    ]
                )
                FRAGMENTS.append(
                    [
                        SINGLE
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - MW("C6H10O5")
                        + 2 * H2O,
                        10,
                        "Glucosyl ceramide - water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        SINGLE
                        - 2 * MW("C11H19NO9")
                        - MW("C8H15NO6")
                        - MW("C6H10O5")
                        + 3 * H2O,
                        22,
                        "Glucosyl ceramide",
                    ]
                )
                FRAGMENTS.append(
                    [
                        SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") + 3 * H2O,
                        59,
                        "Lactosyl ceramide",
                    ]
                )

        return FRAGMENTS


# x  = AcGD2("AcGD2",[[18,1,1], [18,0,0]],  adduct="[M+Na]+")
# print x.printNist()
