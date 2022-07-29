from lipidRules import *

class LysoHexCer(LysoSphingoLipid):   #e.g. psychosine



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]
			FRAGMENTS.append([PREC, 1000, "precursor"])


			if adduct in NEG_ADDUCTS:  
				if adduct == "[M+Cl]-" or adduct == "[M+FA]-" or adduct == "[M+37Cl]-":
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON, 1000, "pre-adduct"] )
				FRAGMENTS.append([MW("C4H6O3")-PROTON, 1000, "hexose fragment"] )
				FRAGMENTS.append([MW("C5H6O3")-PROTON, 1000, "hexose fragment"] )
				FRAGMENTS.append( [MW("C6H10O5")-PROTON, 1000, "deoxyhexose"])
		
			
			return(FRAGMENTS)    




# x  = LysoHexCer("LysoHexCer",[[16,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("LysoHexCer", [[16,0,0]], "[M+H]+")