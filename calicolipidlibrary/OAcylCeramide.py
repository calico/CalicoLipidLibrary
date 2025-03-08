from lipidRules import *

#For the purposes of this class, I am defining the O-Acyl as the sn3 = self.chains[2]
class OAcylCeramide(AcylSphingoLipid):
    pos_adduct_set = ["[M+H]+"]
    def theoreticalDigest(self):
        FRAGMENTS = []
        adduct = self.adduct

        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        if adduct in POS_ADDUCTS:
            if adduct == "[M+H]+":
                FRAGMENTS.append([PREC, 2, "precursor"])
                FRAGMENTS.append([PREC - H2O, 900, "pre-H2O"]) #a
                FRAGMENTS.append([PREC - NL(self.chains[2]), 75, "pre-acyl"]) #b
                FRAGMENTS.append(
                    [PREC - NL(self.chains[2]) - H2O, 750, "pre-acyl-water"] #c
                )
                FRAGMENTS.append([PREC - NL(self.chains[1]) - NL(self.chains[2]), 1000, "sphingoid base"]) #e
                FRAGMENTS.append([PREC - NL(self.chains[1]) - NL(self.chains[2]) - MW("NH3"), 50, "sphingoid base - NH3"])

                FRAGMENTS.append([NL(self.chains[1]) - H2O + MW("NC2H3") + ADDUCT[adduct], 80, "sn2 amide + CH2"]) #d
                FRAGMENTS.append([NL(self.chains[2]) -H2O + MW("NH") + ADDUCT[adduct], 50, "sn2 amide"]) #f

        return FRAGMENTS
