from lipidRules import *

class LPC(LysoGPL):
    neg_adduct_set = GPL.neg_adduct_set + ["[M-CH3]-"]


    def theoreticalDigest(self):
            FRAGMENTS = []
            adduct = self.adduct


            MASS = MW_list(self.MF())
            PREC  = MASS + ADDUCT[self.adduct]
            if adduct in POS_ADDUCTS:

                if adduct == "[M+H]+":  
                    FRAGMENTS.append( [PREC, 295, "precursor"] )
                    FRAGMENTS.append( [MW("C5H14NO"), 758, "Choline"] )
                    FRAGMENTS.append( [MW("C5H14NO")-H2O, 164, "Dehydro Choline"] )                  
                    FRAGMENTS.append( [PREC-H2O, 73, "NL water"] )
                    FRAGMENTS.append( [PREC-MW("C3H9N")-H2O, 9, "NL trimethylamine + water"] )
                    FRAGMENTS.append( [PREC-MW("C5H14NO4P"), 28, "NL phospho choline"] )
                    FRAGMENTS.append( [MW("C8H20NO6P")+PROTON, 13, "glycerophosphocholine"] )        	
                    FRAGMENTS.append( [MW("C5H14NO4P")+ADDUCT[adduct], 1000, "phospocholine"] )
                    FRAGMENTS.append( [MW("C2H6O4P"), 56, "C2H6O4P"] )
                    FRAGMENTS.append( [MW("H3O4P")+ADDUCT[adduct], 2, "phosphoric acid "] )

                if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
                    FRAGMENTS.append( [PREC, 239, "precursor"] )
                    FRAGMENTS.append( [MW("C5H14NO"), 1000, "Choline"] )
                    FRAGMENTS.append( [MW("C5H14NO")-H2O, 195, "Dehydro Choline"] )                             
                    FRAGMENTS.append( [PREC - MW("C3H9N") ,539, "NL trimethylamine"] )
                    FRAGMENTS.append( [PREC-MW("C5H14NO")+PROTON, 13, "NL choline"] )        
                    FRAGMENTS.append( [MW("C3H9O6P")-H2O+ADDUCT[adduct], 9, "glycerol-P - water + adduct"] )
                    FRAGMENTS.append( [MW("C2H5O4P")+ADDUCT[adduct], 243, "C2H5O4P + adduct"] )
                    FRAGMENTS.append( [PREC - MW("C5H14NO4P"), 5, "NL PhosphoCholine"] )
                    FRAGMENTS.append( [PREC - MW("C5H14NO4P") - ADDUCT[adduct] + PROTON, 123, "NL PhosphoCholine + adduct"] )

                
            elif adduct in NEG_ADDUCTS:
                #FRAGMENTS.append([PREC, 295, "precursor"])


                if self.chains[0][0] != 0:
                    i = 0

                    FRAGMENTS.append([PREC - ADDUCT[adduct] - MW("CH3"), 307, "pre-adduct-CH3"])


                    FRAGMENTS.append([MW("HPO3") - PROTON, 35, "PO3-"])
                    FRAGMENTS.append([MW("C5H15NO4P") - PROTON - MW("CH3"), 21, "PC - CH3"])
                    FRAGMENTS.append([MW("C3H9O6P") - PROTON - H2O, 8, "deoxyglycerol phosphate"])

                    FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn1"])
                    FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("CH3") - NL(self.chains[i])+H2O, 7, "NL sn1 ketene"])
                    FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("CH3") - NL(self.chains[i]), 88, "NL sn1"])

                elif self.chains[1][0] != 0:
                    i = 1
                    FRAGMENTS.append([PREC - ADDUCT[adduct] - MW("CH3"), 156, "pre-adduct-CH3"])

                    FRAGMENTS.append([MW("HPO3") - PROTON, 11, "PO3-"])
                    FRAGMENTS.append([MW("C5H15NO4P") - PROTON - MW("CH3"), 5, "PC - CH3"])
                    FRAGMENTS.append([MW("C3H9O6P") - PROTON - H2O, 8, "deoxyglycerol phosphate"])


                    FRAGMENTS.append([NL(self.chains[i]) - PROTON, 1000, "sn1"])
                    FRAGMENTS.append([PREC - ADDUCT[adduct] - MW("CH3") - NL(self.chains[i]) + H2O, 22, "NL sn1 ketene"])
                    FRAGMENTS.append([PREC - ADDUCT[adduct] - MW("CH3") - NL(self.chains[i]), 5, "NL sn1"])



            return(FRAGMENTS)

    @classmethod
    def generateLibrary(cls,target=None, mode="pos"):
        if target: handle = open(target, 'a+')
        if mode == "pos":
            adduct_set = cls.pos_adduct_set
        elif mode =="neg":
            adduct_set = cls.neg_adduct_set
        parent = cls.__bases__[0]
        class_name = cls.__name__
        for c in parent.chain_sets:
            for adduct in adduct_set:
                x  = cls(class_name, c, adduct=adduct)
                content = x.printNist()
                if target: handle.write(content)
                else: sys.stdout.write(content)
        if target: handle.close()





x  = LPC("LPC",[[18,0,0],[0,0,0]],  adduct="[M+H]+")
print x.printNist()

x  = LPC("LPC",[[18,0,0],[0,0,0]],  adduct="[M+Na]+")
print x.printNist()
