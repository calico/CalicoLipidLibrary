from lipidRules import *

class MG(MAG):
		neg_adduct_set = ["[M-H]-"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 134, "precursor"])
					FRAGMENTS.append( [PREC-H2O, 1000, "NL water"] )
					FRAGMENTS.append( [MW("C3H8O3")-H2O+PROTON, 37, "deoxyglycerol"])

					for i in range(0,len(self.chains)):
						chain = str(i+1)

						FRAGMENTS.append( [NL(self.chains[i])+PROTON, 176, "sn"+chain])
						FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 296, "sn"+chain+" acylium"])
						FRAGMENTS.append( [NL(self.chains[i])-2*H2O+PROTON, 53, "sn"+chain+" acylium - water"])

				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:

					FRAGMENTS.append( [MW("C3H8O3")-H2O+ADDUCT[adduct], 37, "deoxyglycerol + adduct"])
					FRAGMENTS.append( [MW("C3H8O3")+ADDUCT[adduct], 37, "glycerol + adduct"])

					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-H2O+PROTON, 1000, "NL adduct + water"] )
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-H2O-H2O+PROTON, 1000, "NL adduct + 2*water"] )
			
					FRAGMENTS.append( [PREC-ADDUCT[adduct]+PROTON, 1000, "NL adduct"] )
					FRAGMENTS.append( [PREC-H2O, 1000, "NL H2O"] )
					FRAGMENTS.append( [PREC-H2O-H2O, 1000, "NL 2*H2O"] )
					for i in range(0,len(self.chains)):
						chain = str(i+1)
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 1000, "sn"+chain+" acylium + adduct"] )
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 1000, "sn"+chain+" acylium"] )
						FRAGMENTS.append( [NL(self.chains[i])+PROTON, 1000, "sn"+chain] )

						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O+MW("C3H6O2"), 1000, "sn"+chain+" + glycerol - water"] )
						FRAGMENTS.append( [NL(self.chains[i])+PROTON+MW("C3H6O2"), 1000, "sn"+chain+" + glycerol"] )

			if adduct in NEG_ADDUCTS:

				if adduct == "[M-H]-":
					FRAGMENTS.append([PREC, 55, "precursor"])
					FRAGMENTS.append([NL(self.chains[0])-PROTON, 1000, "sn1"])


			return(FRAGMENTS)    




# x  = MG("MG",[[14,0,0]],  adduct="[M-H]-")
# print x.printNist()
#
# x  = MG("MG",[[16,0,0]],  adduct="[M-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("MG", [[14,0,0]], "[M-H]-")