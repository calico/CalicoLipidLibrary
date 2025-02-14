from lipidRules import *


class GcGM2(SphingoLipid):
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]
        FRAGMENTS.append([PREC, 1000, "precursor"])

        if adduct in POS_ADDUCTS:
            if adduct == "[M+H]+":
                FRAGMENTS.append([PREC - H2O, 1000, "pre - water"])

                FRAGMENTS.append(
                    [
                        MW("C11H19NO10") - H2O + PROTON,
                        1000,
                        "N-Gc-Neuraminic Acid - water + H+",
                    ]
                )
                FRAGMENTS.append(
                    [
                        MW("C8H15NO6") - H2O + PROTON,
                        1000,
                        "N-Ac-glucosamine - water + H+",
                    ]
                )
                FRAGMENTS.append(
                    [
                        MW("C8H15NO6") - 2 * H2O + PROTON,
                        1000,
                        "N-Ac-glucosamine - 2*water + H+",
                    ]
                )

                FRAGMENTS.append(
                    [PREC - MW("C11H19NO10"), 1000, "pre - N-Gc-Neuraminic Acid"]
                )
                FRAGMENTS.append(
                    [PREC - MW("C8H15NO6"), 1000, "pre - N-Ac-Glucosamine"]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C11H19NO10") - MW("C8H15NO6") + H2O,
                        1000,
                        "pre - N-Gc-Neuraminic Acid - N-Ac-Glucosamine ",
                    ]
                )

                FRAGMENTS.append(
                    [
                        PREC - MW("C11H19NO10") - MW("C8H15NO6") - MW("C6H10O5") + H2O,
                        1000,
                        "pre - N-Gc-Neuraminic Acid - N-Ac-Glucosamine - glucose",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        + H2O,
                        1000,
                        "pre - N-Gc-Neuraminic Acid - N-Ac-Glucosamine - glucose - galactose",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C11H19NO10") - MW("C8H15NO6") - 2 * MW("C6H10O5"),
                        1000,
                        "pre - N-Gc-Neuraminic Acid - N-Ac-Glucosamine - glucose - galactose - water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        - NL(self.chains[1])
                        + H2O * 2,
                        1000,
                        "pre- (head + sn2)",
                    ]
                )

                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        - NL(self.chains[1])
                        + H2O,
                        1000,
                        "pre - (head + sn2 -water )",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        - NL(self.chains[1])
                        + H2O
                        - MW("C"),
                        1000,
                        "pre - (head + sn2 -CH2O )",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[1]) - H2O + MW("NH3") + PROTON,
                        1000,
                        "sn2 fatty amide",
                    ]
                )

            if adduct == "[M+Na]+" or adduct == "[M+K]+":
                FRAGMENTS.append([PREC - MW("CO2"), 1000, "Precursor - CO2"])
                # FRAGMENTS.append( [MW("C11H19NO10")+ADDUCT[adduct], 1000, "pre - N-Gc-Neuraminic Acid"] )
                FRAGMENTS.append(
                    [PREC - MW("C11H19NO10") + H2O, 1000, "pre - N-Gc-Neuraminic Acid"]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C11H19NO10") - MW("C8H15NO6") + 2 * H2O,
                        1000,
                        "pre - N-Gc-Neuraminic Acid - N-Ac-Glucosamine",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - MW("C6H10O5")
                        + 2 * H2O,
                        1000,
                        "pre - N-Gc-Neuraminic Acid - N-Ac-Glucosamine - glucose",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        + 2 * H2O,
                        1000,
                        "pre - N-Gc-Neuraminic Acid - N-Ac-Glucosamine - glucose - galactose",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        + H2O,
                        1000,
                        "pre - N-Gc-Neuraminic Acid - N-Ac-Glucosamine - glucose - galactose - water",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - 2 * MW("C6H10O5")
                        - NL(self.chains[1])
                        + H2O,
                        1000,
                        "pre- (head + sn2)",
                    ]
                )
                # FRAGMENTS.append( [NL(self.chains[1])-H2O+PROTON, 1000, "sn2 acylium ion"])
                # FRAGMENTS.append( [NL(self.chains[1])-H2O-H2O+PROTON, 1000, "sn2 acylium ion - H2O"])

        elif adduct in NEG_ADDUCTS:
            if adduct == "[M-H]-":  # all by analogy to AcGM3
                FRAGMENTS.append([PREC, 1000, "Precursor"])
                FRAGMENTS.append([PREC - H2O, 3, "NL water"])
                FRAGMENTS.append([PREC - MW("CO2"), 5, "NL CO2"])
                FRAGMENTS.append(
                    [MW("C11H19NO10") - PROTON, 15, "N-Ac-Neuraminic Acid"]
                )
                FRAGMENTS.append(
                    [MW("C11H19NO10") - H2O - PROTON, 388, "deoxy-N-Ac-Neuraminic Acid"]
                )

                FRAGMENTS.append(
                    [PREC - MW("C11H19NO10") + H2O, 76, "NL N-Ac-Neuraminic Acid"]
                )
                FRAGMENTS.append(
                    [
                        PREC - MW("C11H19NO10") - MW("C8H15NO6") + 2 * H2O,
                        28,
                        "pre - N-Ac-Neuraminic Acid - N Ac Glucosamine",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - MW("C6H10O5")
                        + 2 * H2O,
                        15,
                        "pre - N-Ac-Neuraminic Acid - N Ac Glucosamine - glucose",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        + 2 * H2O,
                        33,
                        "pre - N-Ac-Neuraminic Acid - N Ac Glucosamine - glucose - galactose",
                    ]
                )
                FRAGMENTS.append(
                    [
                        PREC
                        - MW("C11H19NO10")
                        - MW("C8H15NO6")
                        - 2 * MW("C6H10O5")
                        + H2O,
                        2,
                        "pre - N-Ac-Neuraminic Acid - N Ac Glucosamine - glucose - galactose - water",
                    ]
                )

        return FRAGMENTS


# x  = GcGM2("GcGM2",[[18,1,1], [24,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("GcGM2", [[18,1,1], [24,0,0]], "[M+H]+")
