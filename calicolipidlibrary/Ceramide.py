from lipidRules import *
# need to fix nomenclature
class Ceramide(SphingoLipid):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:



				if adduct == "[M+H]+":

					if self.chains[0][1] == 0: #for sphinganines

						FRAGMENTS.append( [PREC - NL(self.chains[1])+H2O, 92, " NL sn2 ketene"]) #unique to sphinganine

						FRAGMENTS.append([PREC, 101, "precursor"])
						FRAGMENTS.append([PREC - H2O, 1000, "pre-H2O"])
						FRAGMENTS.append( [PREC-H2O-H2O, 65, "pre-H2O-H2O"] )
						FRAGMENTS.append( [PREC - NL(self.chains[1]), 448, "pre-sn2"])
						FRAGMENTS.append( [PREC - NL(self.chains[1])-H2O, 278, "pre-sn2-H2O"])
						FRAGMENTS.append( [PREC - NL(self.chains[1])-MW("CH2O"), 37, "pre-sn2-CH2O"])
						FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 308, "fatty amide"])
						FRAGMENTS.append( [NL(self.chains[1])-H2O+PROTON, 3, "sn2 acylium ion"])


					else:  #for sphingosines
						FRAGMENTS.append([PREC, 2, "precursor"])
						FRAGMENTS.append([PREC - H2O, 101, "pre-H2O"])
						FRAGMENTS.append([PREC - MW("CH2O") - H2O, 13, "pre-CH2O-H2O"])
						FRAGMENTS.append( [PREC-H2O-H2O, 31, "pre-H2O-H2O"] )
						FRAGMENTS.append( [PREC - NL(self.chains[1]), 127, "pre-sn2"])
						FRAGMENTS.append( [PREC - NL(self.chains[1])-H2O, 1000, "pre-sn2-H2O"])
						FRAGMENTS.append( [PREC - NL(self.chains[1])-MW("CH2O"), 101, "pre-sn2-CH2O"])
						FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 11, "fatty amide"])
						FRAGMENTS.append( [NL(self.chains[1])-H2O+PROTON, 1, "sn2 acylium ion"])



				if adduct in ["[M+Li]+", "[M+Na]+", "[M+K]+",  "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 1000, "precursor"])

					FRAGMENTS.append([PREC - H2O, 1000, "pre-H2O"])
					FRAGMENTS.append([PREC - MW("CH2O") - H2O, 1000, "pre-CH2O-H2O"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH")+ADDUCT[adduct], 1000, "sn2 fatty amide + adduct"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O+PROTON, 1000, "sn2 acylium ion"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O-H2O+PROTON, 1000, "sn2 acylium ion - H2O"])
					FRAGMENTS.append( [PREC - NL(self.chains[1])+H2O, 1000, "pre-sn2+H2O+ADDUCT"])

			
			elif adduct in NEG_ADDUCTS:

				if adduct == "[M-H]-":
					FRAGMENTS.append([PREC, 1000, "precursor"])

				elif adduct == "[M+FA]-" or adduct == "[M+Cl]-" or adduct == "[M+37Cl]-":
					FRAGMENTS.append([PREC, 30, "precursor"])
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON, 1000, "NL adduct"] )

				FRAGMENTS.append( [PREC - ADDUCT[adduct] - PROTON - H2O, 20, "pre-H2O"])

				FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - PROTON, 90, "LCB M-H"])
				FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("CH2O") - PROTON, 20, "LCB -CO2H2"])
				FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("CH4O") - PROTON, 62, "LCB -CO2H4"])
				FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("H5NO") - PROTON, 64, "LCB - H2O - NH3"])
				FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - PROTON - MW("CH4O2"), 18, "LCB-COH2-H2O"])
		
				FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("C2H7ON") - PROTON, 13, "N-acyl-ethanolamine - H"])
				FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("C2H5ON") - PROTON, 117, "N-acyl-ethanolamine - H2 - H"])
				FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON, 639, "N-acyl-ethanolamine - H2O - H"])
				FRAGMENTS.append( [ NL(self.chains[1]) - PROTON, 120, "Fatty Acid - H"])
				FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("NH3") - PROTON, 105, "Fatty Amide - H"])
				FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 219, "Fatty Amide - NH3 - H"])

				if (self.chains[1][2] == 1):
					FRAGMENTS.append( [ NL(self.chains[1]) - H2O -  MW("CO") - PROTON, 300, "Fatty Acyl - CH2O2"])
					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("H2")- PROTON, 1000, "Fatty Acyl +H2"])
					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON + MW("H2") - MW("O"), 1000, "Fatty Amide - O + 2H - NH3 - H"])

				if (self.chains[0][2] > 0): 
					FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("H2O")  - MW("NC2H3") - PROTON, 25, "LCB - NC2H3 "])
					FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("H2O") - MW("C2H5N") - PROTON, 190, "LCB-CO2H4-CNH3"])

				if (self.chains[0][1] == 0):
					FRAGMENTS.append([PREC - ADDUCT[adduct] - PROTON - MW("CH4O2"), 223, "pre-COH2-H2O"])
					FRAGMENTS.append([PREC - ADDUCT[adduct] - PROTON - MW("CH4O"), 200, "pre-COH4"])

					FRAGMENTS.append( [PREC - ADDUCT[adduct] - PROTON - MW("CH2O"), 2, "pre-COH2"])

				if (self.chains[0][1] == 1):
					FRAGMENTS.append([PREC - ADDUCT[adduct] - PROTON - MW("CH4O2"), 74, "pre-COH2-H2O"])
					FRAGMENTS.append([PREC - ADDUCT[adduct] - PROTON - MW("CH4O"), 67, "pre-COH4"])

					FRAGMENTS.append( [PREC - ADDUCT[adduct] - PROTON - MW("CH2O"), 110, "pre-COH2"])

				if (self.chains[0][1] == 2):
					FRAGMENTS.append([PREC - ADDUCT[adduct] - PROTON - MW("CH4O2"), 105, "pre-COH2-H2O"])
					FRAGMENTS.append([PREC - ADDUCT[adduct] - PROTON - MW("CH4O"), 19, "pre-COH4"])

				if (self.chains[0][2] == 2):
					FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("H2O") - MW("C3H7NO") - PROTON, 50, "LCB-CO2H4-CNH3 - (CH2O if t) "])  #50 is guess, not from spectrum
					FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("H2O")  - PROTON, 55, "LCB - water "])


				
			
			
			return(FRAGMENTS)    

# calicolipidlibrary.print_spectrum("Ceramide",[[18,1,1], [16,0,0]], "[M+H]+")

# x  = Ceramide("Ceramide",[[18,1,1], [16,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# x  = Ceramide("Ceramide",[[18,1,1], [24,1,0]],  adduct="[M+H]+")
# print x.printNist()
#
# x  = Ceramide("Ceramide",[[18,0,1], [15,0,0]],  adduct="[M-H]-")
# print x.printNist()
