from lipidRules import *

class Ceramide_P(SphingoLipid):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]


				
			if adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])

				FRAGMENTS.append([MW("HPO3")-PROTON, 800, "phosphite"] )
				FRAGMENTS.append([MW("H3PO4")-PROTON, 230, "phosphate"] )
				FRAGMENTS.append([PREC - NL(self.chains[1]) + H2O, 7, "NL acyl ketene"] )
				FRAGMENTS.append([PREC - NL(self.chains[1]), 2, "NL acyl ketene + water"] )


  			if adduct in POS_ADDUCTS:  
				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 1000, "precursor"])

					FRAGMENTS.append([MW("H3PO4")+ADDUCT[adduct], 95, "Phosphoric acid + adduct"] )
					FRAGMENTS.append([PREC - MW("H3PO4"), 995, "NL Phosphoric acid"] )					
					FRAGMENTS.append([PREC - MW("H3PO4")-ADDUCT[adduct]+PROTON-H2O, 76, "NL Phosphoric acid + water + adduct"] )					

					FRAGMENTS.append([NL(self.chains[1])-H2O + MW("NH3")+PROTON, 95, "fatty amide"] )		
					FRAGMENTS.append( [PREC - NL(self.chains[1]) - MW("H3O4P")  - ADDUCT[adduct] + PROTON, 7, "N'': NL PhosphoCholine + sn2 ketene + water + ADDUCT"] )  

					FRAGMENTS.append([NL(self.chains[1])-H2O + MW("NH3")+MW("C2") +PROTON, 44, "fatty amide + C2"] )		


				if adduct == "[M+H]+":
			
					FRAGMENTS.append([PREC, 1, "precursor"])
					FRAGMENTS.append([PREC, 224, "NL water"])
					FRAGMENTS.append([PREC - MW("H3PO4"), 136, "NL Phosphoric acid"] )					
					FRAGMENTS.append([PREC - MW("H3PO4")-ADDUCT[adduct]+PROTON-H2O, 212, "NL Phosphoric acid + water + adduct"] )					

					FRAGMENTS.append([NL(self.chains[1])-H2O + MW("NH3")+PROTON, 26, "fatty amide"] )		
					FRAGMENTS.append( [PREC - NL(self.chains[1]) - MW("H3O4P")  - ADDUCT[adduct] + PROTON, 1000, "N'': NL PhosphoCholine + acyl ketene + water + ADDUCT"] )  
					FRAGMENTS.append( [PREC - NL(self.chains[1]) - MW("H3O4P") - MW("NH3")  - ADDUCT[adduct] + PROTON, 8, "NL PhosphoCholine + fatty amide + water + ADDUCT"] )  



					FRAGMENTS.append([NL(self.chains[1])-H2O + MW("NH3")+MW("C2") +PROTON, 33, "fatty amide +C2"] )		



			
			return(FRAGMENTS)    




# x  = Ceramide_P("Ceramide_P",[[18,1,1], [12,0,0]],  adduct="[M-H]-")
# print x.printNist()
#
#
# x  = Ceramide_P("Ceramide_P",[[18,1,1], [24,0,0]],  adduct="[M-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("Ceramide_P", [[18,1,1], [12,0,0]], "[M-H]-")