from lipidRules import *

class PC(GPL):

	neg_adduct_set = GPL.neg_adduct_set + ["[M-CH3]-"]

	def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
			
				if adduct == "[M+H]+": 
					FRAGMENTS.append([PREC, 990, "precursor"])

					FRAGMENTS.append( [MW("C2H5O4P")+PROTON, 15, "C2H5O4P+Proton"] ) 
					FRAGMENTS.append( [MW("C5H14NO4P")+ADDUCT[adduct], 649, "Phospocholine"] )      	
					FRAGMENTS.append( [PREC - MW("C3H9N"), 830, "NL Trimethylamine"] )       
			
			
					FRAGMENTS.append( [PREC-MW("C5H15NO4P")+PROTON, 1000, "NL choline-P"] )
					FRAGMENTS.append( [MW("C5H13NO")+ADDUCT[adduct], 5, "Choline"] )
					FRAGMENTS.append( [MW("H3O4P")+ADDUCT[adduct], 1, "PO4H4+"] )

				
					i = 0
					chain = str(i+1)
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - ADDUCT[adduct] + PROTON, 7, "NL sn" + chain + " -ADDUCT + PROTON"] )
					FRAGMENTS.append( [PREC - NL(self.chains[i]) , 7, "NL sn" + chain ] )
					FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 10, "sn" + chain + " acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - MW("C3H9N"), 33, "NL sn" + chain + "-trimethylamine"])
	
		
					i = 1
					chain = str(i+1)
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - ADDUCT[adduct] + PROTON, 7, "NL sn" + chain + " -ADDUCT + PROTON"] )
					FRAGMENTS.append( [PREC - NL(self.chains[i]) , 2, "NL sn" + chain ] )
					FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 13, "sn" + chain + " acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - MW("C3H9N"), 12, "NL sn" + chain + "-trimethylamine"])


				
				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 1000, "precursor"])

					FRAGMENTS.append( [MW("C2H5O4P")+PROTON, 2, "C2H5O4P+Proton"] ) 

					FRAGMENTS.append( [MW("C2H5O4P")+ADDUCT[adduct], 567, "Phos + C2 + adduct"] )
					FRAGMENTS.append( [PREC - MW("C3H9N"), 409, "NL Trimethylamine"] )       
					FRAGMENTS.append( [PREC - (MW("C5H15NO4P")-PROTON), 464, "NL PhosphoCholine"] )
					FRAGMENTS.append( [PREC - (MW("C5H15NO4P") + ADDUCT[adduct])+2*PROTON, 207, "NL (PhosphoCholine + ADDUCT) +2*PROTON"] ) #should 2x proton really be 2x Hydrogen?
					FRAGMENTS.append( [(MW("C5H14NO4P")+PROTON), 32, "PhosphoCholine + PROTON"] )       
		
					i = 0
					chain = str(i+1)
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - ADDUCT[adduct] + PROTON, 7, "NL sn" + chain + " -ADDUCT + PROTON"] )
					FRAGMENTS.append( [PREC - NL(self.chains[i]) , 7, "NL sn" + chain ] )
					FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 10, "sn" + chain + " acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - MW("C3H9N"), 33, "NL sn" + chain + "trimethylamine"])
	
		
					i = 1
					chain = str(i+1)
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - ADDUCT[adduct] + PROTON, 7, "NL sn" + chain + " -ADDUCT + PROTON"] )
					FRAGMENTS.append( [PREC - NL(self.chains[i]) , 2, "NL sn" + chain ] )
					FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 13, "sn" + chain + " acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - MW("C3H9N"), 12, "NL sn" + chain + "trimethylamine"])




		
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1, "precursor"])
				if adduct != "[M-CH3]-":
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("CH3"), 255, "NL adduct  + CH3"] )
				#FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C5H12N"), 1000, "pre-adduct-choline"] )   	#choline = C5H14NO+  MAss is off by ~0.15 D
				#FRAGMENTS.append( [MW("C3H9O6P")-H2O+ADDUCT[adduct], 1000,  "glycerol3P-H2O"] )
				#FRAGMENTS.append( [MW("H3PO4")-PROTON,  1000,  "H2PO4-"] )
				FRAGMENTS.append( [MW("HPO3")-PROTON,  26,  "PO3-"] )
				FRAGMENTS.append( [MW("C5H15NO4P")-PROTON-MW("CH3"),  35,  "PC - CH3"] )
				FRAGMENTS.append( [MW("C8H21NO6P")-PROTON-H2O-MW("CH3"),  26,  "gPC- CH3-H2O"] )
				i = 0
				chain = str(i+1)
				FRAGMENTS.append( [NL(self.chains[i])-PROTON, 136, "sn"+chain])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("CH3") - NL(self.chains[i])+H2O, 10, "NL sn"+chain + "+H2O"])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("CH3") - NL(self.chains[i]), 7, "NL sn"+chain])

				i = 1
				chain = str(i+1)
				FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn"+chain])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("CH3") - NL(self.chains[i])+H2O, 25, "NL sn"+chain + "+H2O"])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("CH3") - NL(self.chains[i]), 5, "NL sn"+chain])

			return(FRAGMENTS)    




# x  = PC("PC",[[16,0,0], [18,1,0]],  adduct="[M+Na]+")
# print x.printNist()
#
#
#
# x  = PC("PC",[[16,0,0], [18,1,0]],  adduct="[M+H]+")
# print x.printNist()
#
#
# x  = PC("PC",[[16,0,0], [16,0,0]],  adduct="[M+FA-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("PC",[[18,0,0], [20,4,0]],  adduct="[M+FA-H]-")