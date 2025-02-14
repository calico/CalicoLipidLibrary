from lipidRules import *


class CL(CardioLipin):
    neg_adduct_set = ["[M-H]-", "[M-2H]2-"]

    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = (MASS + ADDUCT[self.adduct]) / abs(ADDUCT_CHARGE[self.adduct])
        SC_PREC = MASS - PROTON

        if adduct in NEG_ADDUCTS:
            if adduct == "[M-H]-":
                FRAGMENTS.append([PREC, 703, "precursor"])

                FRAGMENTS.append([MW("C3H9O6P") - PROTON, 9, "Glycerol-P"])
                FRAGMENTS.append(
                    [MW("C3H9O6P") - H2O - PROTON, 463, "Glycerol-P - water"]
                )
                for i in [0, 2]:
                    chain = str(i + 1)

                    FRAGMENTS.append([NL(self.chains[i]) - PROTON, 500, "sn" + chain])
                    FRAGMENTS.append([PREC - NL(self.chains[i]), 1, "nL sn" + chain])
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + MW("C3H7O5P") - PROTON,
                            42,
                            "sn" + chain + " LPA",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O + MW("C3H7O5P") - PROTON,
                            99,
                            "sn" + chain + " LPA - water",
                        ]
                    )

                for i in [1, 3]:
                    chain = str(i + 1)

                    FRAGMENTS.append([NL(self.chains[i]) - PROTON, 325, "sn" + chain])
                    FRAGMENTS.append([PREC - NL(self.chains[i]), 2, "nL sn" + chain])
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + MW("C3H7O5P") - PROTON,
                            14,
                            "sn" + chain + " LPA",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O + MW("C3H7O5P") - PROTON,
                            41,
                            "sn" + chain + " LPA - water",
                        ]
                    )

                FRAGMENTS.append(
                    [
                        NL(self.chains[0])
                        + NL(self.chains[1])
                        + MW("C3H5O4P")
                        - PROTON,
                        175,
                        "PA: sn1+sn2",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[2])
                        + NL(self.chains[3])
                        + MW("C3H5O4P")
                        - PROTON,
                        175,
                        "PA: sn3+sn4",
                    ]
                )

                FRAGMENTS.append(
                    [
                        NL(self.chains[0])
                        + NL(self.chains[1])
                        + MW("C3H5O4P")
                        + MW("C3H6O2")
                        - H2O
                        - PROTON,
                        15,
                        "PG: sn1+sn2",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[2])
                        + NL(self.chains[3])
                        + MW("C3H5O4P")
                        + MW("C3H6O2")
                        - H2O
                        - PROTON,
                        15,
                        "PG: sn3+sn4",
                    ]
                )

                FRAGMENTS.append(
                    [
                        NL(self.chains[0])
                        + NL(self.chains[1])
                        + MW("C3H5O4P")
                        + MW("C3H6O2")
                        + MW("PO3H")
                        - H2O
                        - PROTON,
                        19,
                        "phosphoPG: sn1+sn2",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[2])
                        + NL(self.chains[3])
                        + MW("C3H5O4P")
                        + MW("C3H6O2")
                        + MW("PO3H")
                        - H2O
                        - PROTON,
                        19,
                        "phosphoPG: sn3+sn4",
                    ]
                )

            if adduct == "[M-2H]2-":
                FRAGMENTS.append([PREC, 63, "precursor"])

                FRAGMENTS.append([MW("C3H9O6P") - PROTON, 1, "Glycerol-P"])
                FRAGMENTS.append(
                    [MW("C3H9O6P") - H2O - PROTON, 149, "Glycerol-P - water"]
                )
                for i in [0, 2]:
                    chain = str(i + 1)

                    FRAGMENTS.append([NL(self.chains[i]) - PROTON, 430, "sn" + chain])
                    FRAGMENTS.append(
                        [SC_PREC - NL(self.chains[i]), 11, "nL sn" + chain]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + MW("C3H7O5P") - PROTON,
                            2,
                            "sn" + chain + " LPA",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O + MW("C3H7O5P") - PROTON,
                            73,
                            "sn" + chain + " LPA - water",
                        ]
                    )

                    FRAGMENTS.append(
                        [
                            PREC - (NL(self.chains[i]) - H2O) / 2,
                            7,
                            "nL sn" + chain + " ketene 2-",
                        ]
                    )

                for i in [1, 3]:
                    chain = str(i + 1)

                    FRAGMENTS.append([NL(self.chains[i]) - PROTON, 500, "sn" + chain])
                    FRAGMENTS.append(
                        [SC_PREC - NL(self.chains[i]), 26, "nL sn" + chain]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) + MW("C3H7O5P") - PROTON,
                            42,
                            "sn" + chain + " LPA",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            NL(self.chains[i]) - H2O + MW("C3H7O5P") - PROTON,
                            73,
                            "sn" + chain + " LPA - water",
                        ]
                    )
                    FRAGMENTS.append(
                        [
                            PREC - (NL(self.chains[i]) - H2O) / 2,
                            18,
                            "nL sn" + chain + " ketene 2-",
                        ]
                    )

                FRAGMENTS.append(
                    [
                        NL(self.chains[0])
                        + NL(self.chains[1])
                        + MW("C3H5O4P")
                        - PROTON,
                        7,
                        "PA: sn1+sn2",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[2])
                        + NL(self.chains[3])
                        + MW("C3H5O4P")
                        - PROTON,
                        7,
                        "PA: sn3+sn4",
                    ]
                )

                FRAGMENTS.append(
                    [
                        NL(self.chains[0])
                        + NL(self.chains[1])
                        + MW("C3H5O4P")
                        + MW("C3H6O2")
                        - H2O
                        - PROTON,
                        10,
                        "PG: sn1+sn2",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[2])
                        + NL(self.chains[3])
                        + MW("C3H5O4P")
                        + MW("C3H6O2")
                        - H2O
                        - PROTON,
                        10,
                        "PG: sn3+sn4",
                    ]
                )

                FRAGMENTS.append(
                    [
                        NL(self.chains[0])
                        + NL(self.chains[1])
                        + MW("C3H5O4P")
                        + MW("C3H6O2")
                        + MW("PO3H")
                        - H2O
                        - PROTON,
                        1,
                        "phosphoPG: sn1+sn2",
                    ]
                )
                FRAGMENTS.append(
                    [
                        NL(self.chains[2])
                        + NL(self.chains[3])
                        + MW("C3H5O4P")
                        + MW("C3H6O2")
                        + MW("PO3H")
                        - H2O
                        - PROTON,
                        1,
                        "phosphoPG: sn3+sn4",
                    ]
                )

        return FRAGMENTS


# x  = CL("CL",[[16,0,0], [18,1,0], [16,0,0], [18,1,0]],  adduct="[M-H]-")
# print x.printNist()
#
# x  = CL("CL",[[16,0,0], [18,1,0], [16,0,0], [18,1,0]],  adduct="[M-2H]2-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("CL", [[16,0,0], [18,1,0], [16,1,0], [18,0,0]], "[M-2H]2-")
