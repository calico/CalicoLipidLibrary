import sys
from lipidRules import *

class TG_DI(Triglyceride):

        pos_adduct_set = ["[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+Li]+"]

        def theoreticalDigest(self):
            FRAGMENTS = []
            adduct = self.adduct

            MASS = MW_list(self.MF())
            PREC  = MASS + ADDUCT[self.adduct]

            if adduct in POS_ADDUCTS:

                if adduct == "[M+H]+":
                    FRAGMENTS.append([PREC, 1000, "precursor"])

                if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
                    FRAGMENTS.append([PREC, 1000, "precursor"])

                    chain_list = range(0,len(self.chains))
                    for i in chain_list:
                        chain = str(i+1)
                        MS3_PREC = str(PREC - NL(self.chains[i])-ADDUCT[adduct]+PROTON)
                        ms3_chain_list = list(chain_list)
                        ms3_chain_list.remove(i)
                        FRAGMENTS.append( [PREC - NL(self.chains[i])-ADDUCT[adduct]+PROTON, 61, "NL sn"+chain+" + adduct"])
                        FRAGMENTS.append( [PREC - NL(self.chains[i]),56 , "NL sn"+chain])
                        FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 10, "sn"+chain+" acylium ion"] )
                        FRAGMENTS.append( [PREC - NL(self.chains[i])-ADDUCT[adduct]+PROTON,158 , "ms3-{"+MS3_PREC+"}-MS3 precursor from NL sn"+chain+" + adduct"])
                        for j in ms3_chain_list:
                            ms3_chain = str(j+1)
                            FRAGMENTS.append([NL(self.chains[j]) + PROTON - H2O, 608, "ms3-{"+MS3_PREC+"} sn" + ms3_chain + " acylium ion"])
                            FRAGMENTS.append([NL(self.chains[j]) + PROTON +MW("C3H6O"), 300, "ms3-{"+MS3_PREC+"} sn" + ms3_chain + " + glycerol"])

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
                    x  = cls("TG", c, adduct=adduct)
                    content = x.printNist()
                    if target: handle.write(content)
                    else: sys.stdout.write(content)
            if target: handle.close()





x  = TG_DI("TG",[[17,1,0], [17,1,0], [17,1,0]],  adduct="[M+NH4]+")
print x.printNist()

x  = TG_DI("TG",[[17,1,0], [17,1,0], [17,1,0]],  adduct="[M+Na]+")
print x.printNist()

x  = TG_DI("TG",[[17,1,0], [17,1,0], [17,1,0]],  adduct="[M+K]+")
print x.printNist()

x  = TG_DI("TG",[[17,1,0], [17,1,0], [17,1,0]],  adduct="[M+Li]+")
print x.printNist()



x  = TG_DI("TG",[[15,0,0], [18,1,0], [15,0,0]],  adduct="[M+NH4]+")
print x.printNist()


x  = TG_DI("TG",[[15,0,0], [18,1,0], [15,0,0]],  adduct="[M+Na]+")
print x.printNist()


x  = TG_DI("TG",[[15,0,0], [18,1,0], [15,0,0]],  adduct="[M+K]+")
print x.printNist()


x  = TG_DI("TG",[[15,0,0], [18,1,0], [15,0,0]],  adduct="[M+Li]+")
print x.printNist()


