from lipidRules import *

class FA(singleAcyl):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]
			FRAGMENTS.append([PREC, 1000, "precursor"])

			if adduct in NEG_ADDUCTS:

				FRAGMENTS.append( [PREC-H2O,1000,"loss of water"])

				if self.chains[0][1] > 2:
					FRAGMENTS.append( [PREC-MW("CO2"),1000,"loss of CO2"])
			
			
				


			return(FRAGMENTS)    




# x  = FA("FA",[[16,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("FA", [[16,0,0]], "[M+H]+")