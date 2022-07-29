from lipidRules import *

class Sulfatide(SphingoLipid):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]



			
			if adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])			
				FRAGMENTS.append( [MW("H2SO4")-PROTON, 25, "sulfate"] )
				FRAGMENTS.append( [MW("C6H10O8S")-PROTON , 5, "hexose_sulfate - water"] )
				FRAGMENTS.append( [MW("C6H12O9S")-PROTON , 1, "hexose_sulfate"] )
			
				FRAGMENTS.append([PREC-NL(self.chains[1])+H2O, 5, "NL  ketene"])			
				FRAGMENTS.append([PREC-NL(self.chains[1]), 3, "NL  ketene + water"])			
				FRAGMENTS.append([PREC-MW("C6H10O8S"), 1000, "precursor - hexose-sulfate"])			
				if (self.chains[1][1] == 0):  #if no hydroxyl on acyl group
					FRAGMENTS.append( [MW("C6H12O9S")+MW("C2H3N")-PROTON , 4, "amino-ethenyl hexose_sulfate"] )

		
			return(FRAGMENTS)    




# x  = Sulfatide("Sulfatide",[[18,1,1], [18,0,1]],  adduct="[M-H]-")
# print x.printNist()
#
#
# x  = Sulfatide("Sulfatide",[[18,0,1], [16,1,0]],  adduct="[M+FA-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("Sulfatide", [[18,0,1], [16,1,0]], "[M+FA-H]-")