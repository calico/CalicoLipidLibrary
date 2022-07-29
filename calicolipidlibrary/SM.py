from lipidRules import *

class SM(SphingoLipid):

	neg_adduct_set = GPL.neg_adduct_set + ["[M-CH3]-"]

	def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
			


		
				if adduct == "[M+H]+":  
					FRAGMENTS.append([PREC, 48, "precursor"])				
					FRAGMENTS.append( [MW("C5H13NO")+PROTON, 8, "Choline"] )
					FRAGMENTS.append( [MW("C5H13NO")+PROTON-H2O, 39, "Dehydro Choline"] )				
					FRAGMENTS.append( [PREC-H2O, 6, "precursor - H2O"] )   	
					FRAGMENTS.append( [MW("C5H14NO4P")+PROTON, 1000, "Phospocholine"] )
					FRAGMENTS.append( [MW("C2H6O4P"), 25, "C2H6O4P"] )
					FRAGMENTS.append( [MW("PO4H3")+ADDUCT[adduct], 1000, "phosphoric acid + adduct"] )
					FRAGMENTS.append( [MW_list([self.chains[0][0],(2*self.chains[0][0]-2*self.chains[0][1]-1),1,self.chains[0][2],0,0])+ADDUCT[adduct], 2, "LCB - water"] )
					if self.chains[0][2] > 0: #if not monohydroxy LCB
						FRAGMENTS.append([MW_list([self.chains[0][0], (2 * self.chains[0][0] - 2 * self.chains[0][1] - 1), 1, self.chains[0][2], 0, 0])-H2O + ADDUCT[adduct], 1, "LCB - 2xwater"])


				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 945, "precursor"])
					FRAGMENTS.append( [MW("C5H13NO")+PROTON, 5, "Choline"] )
					FRAGMENTS.append( [MW("C5H13NO")+PROTON-H2O, 181, "Dehydro Choline"] )
					FRAGMENTS.append( [MW("C5H14NO4P")+ADDUCT[adduct], 2, "Phospocholine + ADDUCT"] )
					FRAGMENTS.append( [MW("C2H5O4P")+ADDUCT[adduct], 250, "C2H5O4P + adduct"] )
					FRAGMENTS.append( [PREC - MW("C3H9N"), 837, "NL trimethylamine"] )       
					FRAGMENTS.append( [PREC - MW("C5H14NO4P"), 902, "NL PhosphoCholine"] )       
					FRAGMENTS.append( [PREC - MW("C5H14NO4P") - (ADDUCT[adduct]-PROTON), 6, "NL PhosphoCholine + adduct"] )
					FRAGMENTS.append( [PREC - MW("C5H14NO4P") - (ADDUCT[adduct]-PROTON+H2O), 29, "NL DeoxyPhosphoCholine+ adduct"] )
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("C2H3N")+PROTON, 14, "Fatty amide + C2"])


		
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1, "precursor"])

				FRAGMENTS.append([PREC - ADDUCT[adduct] - MW("CH3"), 1000, "pre-adduct-CH3"])
				FRAGMENTS.append([MW("HPO3") - PROTON, 510, "PO3-"])
				FRAGMENTS.append([MW("C5H15NO4P") - PROTON - MW("CH3"), 507, "phosphocholine - CH3"])
				FRAGMENTS.append([MW("C2H5PO4") - PROTON, 14, "C2H4 + phosphate"])



			return(FRAGMENTS)    



# x  = SM("SM",[[18,1,0], [20,0,0]],  adduct="[M+Na]+")
# print x.printNist()
# x  = SM("SM",[[18,1,1], [20,0,0]],  adduct="[M+Na]+")
# print x.printNist()
# x  = SM("SM",[[18,1,2], [20,0,0]],  adduct="[M+Na]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("SM",[[18,1,0], [20,4,3]],  adduct="[M+Na]+")