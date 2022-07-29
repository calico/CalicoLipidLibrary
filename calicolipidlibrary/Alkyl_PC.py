from lipidRules import *

class Alkyl_PC(alkylGPL):
	neg_adduct_set = GPL.neg_adduct_set + ["[M-CH3]-"]
	def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
			


				if adduct == "[M+H]+":  
					FRAGMENTS.append([PREC, 180, "precursor"])
				
					FRAGMENTS.append( [MW("C5H14NO4P")+ADDUCT[adduct], 1000, "Phospocholine"] )      	
					FRAGMENTS.append( [MW("C2H5O4P")+PROTON, 25, "C2H6O4P"] )
					FRAGMENTS.append( [MW("C5H13NO")+ADDUCT[adduct], 5, "Choline"] )
					FRAGMENTS.append( [MW("H3O4P")+ADDUCT[adduct], 2, "phosphoric acid + adduct"] )
					FRAGMENTS.append( [PREC - NL(self.chains[1])+H2O, 6, "NL sn2 + water"])
	
				
				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 1000, "precursor"])

					FRAGMENTS.append( [MW("C2H5O4P")+ADDUCT[adduct], 432, "C2H5O4P + adduct"] )
					FRAGMENTS.append( [PREC - MW("C3H9N"), 695, "NL Trimethylamine"] )       
					FRAGMENTS.append( [PREC-MW("C5H14NO4P"), 86, "NL phosphocholine"] )
					FRAGMENTS.append( [PREC - MW("C5H14NO4P") - ADDUCT[adduct] + PROTON, 69, "NL PhosphoCholine + ADDUCT"] )
					
					FRAGMENTS.append( [NL(self.chains[0])-H2O-H2O+PROTON, 5, "sn1 ether-linked hydrocarbon"])  
					FRAGMENTS.append( [PREC - NL(self.chains[1]) - ADDUCT[adduct] + PROTON, 6, "NL sn2 + adduct"] )
					FRAGMENTS.append( [PREC - NL(self.chains[1]) - ADDUCT[adduct] + PROTON - MW("C5H14NO4P"), 144, "NL sn2 + phosphocholine + adduct"] )
					FRAGMENTS.append( [NL(self.chains[1])-H2O+PROTON, 7, "sn2 acylium ion"])

			
			elif adduct in NEG_ADDUCTS: 
				if adduct != "[M-H]-":
					FRAGMENTS.append([PREC, 3, "precursor"])

					FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("CH3"), 310, "NL adduct + CH3"] )

					#FRAGMENTS.append( [MW("H3PO4")-PROTON,  1000,  "H2PO4-"] )
					FRAGMENTS.append( [MW("HPO3")-PROTON,  18,  "phosphite"] )
					FRAGMENTS.append( [MW("C5H15NO4P")-PROTON-MW("CH3"),  1,  "PC - CH3"] )

					FRAGMENTS.append( [NL(self.chains[1])-PROTON, 1000, "sn2"])
					#FRAGMENTS.append( [NL(self.chains[1])-H2O-PROTON, 1000, "sn2 ketene"])
					FRAGMENTS.append( [NL(self.chains[0])-MW("O2")+MW("H4")-PROTON, 57, "sn1 alkene"])

					FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C5H12N") - NL(self.chains[1])+H2O , 1, "NL deoxycholine + sn2 + H2O"])
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C5H12N") - NL(self.chains[1]), 14, "NL deoxycholine + sn2"])

					FRAGMENTS.append([PREC - ADDUCT[adduct] - MW("CH3") - NL(self.chains[1]) + H2O, 57, "NL adduct + CH3 + sn2 + H2O"])
					FRAGMENTS.append([PREC - ADDUCT[adduct] - MW("CH3") - NL(self.chains[1]), 25, "NL adduct + CH3 + sn2"])


			return(FRAGMENTS)    




# x  = Alkyl_PC("Alkyl_PC",[[16,0,0], [16,0,0]],  adduct="[M+FA-H]-")
# print x.printNist()
#
#
# x  = Alkyl_PC("Alkyl_PC",[[16,1,0], [18,0,0]],  adduct="[M+FA-H]-")
# print x.printNist()

# calicolipidlibrary.print_spectrum("Alkyl_PC", [[16,1,0], [18,0,0]], "[M+FA-H]-")