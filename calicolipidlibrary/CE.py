from lipidRules import *

class CE(singleAcyl):

		pos_adduct_set = ["[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+Li]+"]

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
					
				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:

					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append( [MW("C27H44") + PROTON,8,"Cholesterol - water"]) 
					FRAGMENTS.append([NL(self.chains[0])+ADDUCT[adduct], 249, "Fatty acyl + adduct"]) 
					FRAGMENTS.append([MW("C11H14")+PROTON, 8, "C11H15+"])
					if self.chains[0][2] > 0:
						FRAGMENTS.append([NL(self.chains[0])+ADDUCT[adduct]-H2O, 87, "Fatty acyl + adduct - H2O"]) 

			return(FRAGMENTS)    




# x  = CE("CE",[[18,2,1]],  adduct="[M+Na]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("CE", [[18,2,1]],"[M+Na]+")
