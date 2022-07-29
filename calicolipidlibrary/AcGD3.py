from lipidRules import *

class AcGD3(SphingoLipid):

		neg_adduct_set = ["[M-H]-", "[M-2H]2-"]
		pos_adduct_set = ["[M+H]+"]

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			MASS = MW_list(self.MF())
			PREC  = (MASS + ADDUCT[self.adduct])/abs(ADDUCT_CHARGE[self.adduct])

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 28, "precursor"])
					FRAGMENTS.append( [PREC-H2O, 14, "NL water"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 9, "NL deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9"), 55, "NL N-Ac-Neuraminic Acid"] )

					FRAGMENTS.append( [PREC - 2*MW("C11H19NO9")+H2O, 155, "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C6H12O6")+2*H2O, 129, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose "] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-2*MW("C6H12O6")+3*H2O, 396, "Ceramide"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-2*MW("C6H12O6")+2*H2O, 104, "Ceramide - water"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-2*MW("C6H12O6")+2*H2O - NL(self.chains[1]) + 2*H2O, 108, "sphingoid base - water"] )

					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-2*MW("C6H12O6")+2*H2O - NL(self.chains[1]) + H2O, 1000, "sphingoid base - water"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-2*MW("C6H12O6")+2*H2O - NL(self.chains[1]) +2*H2O - MW("CH2O"), 90, "sphingoid base - CH2O"] )

					FRAGMENTS.append([MW("C11H19NO9") - H2O + PROTON, 416, "deoxy-N-Ac-Neuraminic Acid"])
					#FRAGMENTS.append([MW("C11H19NO9") + PROTON, 55, "N-Ac-Neuraminic Acid"])
					FRAGMENTS.append( [MW("C11H19NO9")-2*H2O+PROTON, 857, "deoxy-N-Ac-Neuraminic Acid - water "] )

					FRAGMENTS.append( [NL(self.chains[1]) -H2O +MW("NH3")+PROTON, 14, "fatty amide"] )



			# omitted [M=2H]2+ because reference spectrum was bad


			elif adduct in NEG_ADDUCTS:  
				if adduct == "[M+Cl]-" or adduct == "[M+FA]-":
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON, 1000, "NL adduct"] )
				if adduct == "[M-H]-":
					FRAGMENTS.append([PREC, 258, "precursor"])
					FRAGMENTS.append( [PREC - H2O, 11, "NL water"] )
					FRAGMENTS.append( [PREC - MW("CO2"), 59, "NL CO2"] )
					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON, 8, "N-acyl-ethanolamine - H2O"])

					FRAGMENTS.append( [MW("C11H19NO9")-PROTON, 7, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O-PROTON, 1000, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-2*H2O-PROTON, 490, "2x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-2*H2O-PROTON, 161, "2x deoxy-N-Ac-Neuraminic Acid - CO2"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-3*H2O-PROTON, 32, "2x deoxy-N-Ac-Neuraminic Acid - CO2 - water"] )

					FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 792, "NL deoxy-N-Ac-Neuraminic Acid"] )
					#FRAGMENTS.append( [PREC-MW("C11H19NO9"), 4, "NL N-Ac-Neuraminic Acid"] )
					#FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("CO2")+H2O, 2, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )
					FRAGMENTS.append( [PREC-2*MW("C11H19NO9")+2*H2O, 3, "NL 2x deoxy-N-Ac-Neuraminic Acid"] )


					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 2, "Fatty Amide - NH3 - H"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - 2 * MW("C6H10O5") + 2 * H2O, 21, "Ceramide"])
					#FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C6H10O5") + 1 * H2O, 2, "Glucosyl ceramide - water"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C6H10O5") + 2* H2O, 14, "Glucosyl ceramide"])

				if adduct == "[M-2H]2-":
					SINGLE = MASS - PROTON
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append([PREC - H2O / 2, 2, "NL water 2-"])
					#FRAGMENTS.append([PREC - MW("CO2") / 2, 14, "NL CO2 2-"])

					#FRAGMENTS.append([NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON, 16, "N-acyl-ethanolamine - H2O"])
					#FRAGMENTS.append([MW("C11H19NO9") - PROTON, 16, "N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9") - H2O - PROTON, 649, "deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9") + MW("C11H19NO9") - 2 * H2O - PROTON, 208, "2x deoxy-N-Ac-Neuraminic Acid"])
					#FRAGMENTS.append([MW("C11H19NO9") + MW("C11H19NO9") - MW("CO2") - 2 * H2O - PROTON, 11, "2x deoxy-N-Ac-Neuraminic Acid - CO2"])
					#FRAGMENTS.append([MW("C11H19NO9") + MW("C11H19NO9") - MW("CO2") - 3 * H2O - PROTON, 4,  "2x deoxy-N-Ac-Neuraminic Acid - CO2 - water"])

					FRAGMENTS.append([SINGLE - MW("C11H19NO9") + H2O, 88, "NL deoxy-N-Ac-Neuraminic Acid"])


					# FRAGMENTS.append( [SINGLE-MW("C11H19NO9"), 1, "NL N-Ac-Neuraminic Acid"] )
					# FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("CO2")+H2O, 2, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") + 2 * H2O, 109, "NL 2x deoxy-N-Ac-Neuraminic Acid"])

					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 2, "Fatty Amide - NH3 - H"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - 2 * MW("C6H10O5") + 2 * H2O, 21, "Ceramide"])
					#FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C6H10O5") + 2* H2O, 14, "Glucosyl ceramide"])


			return(FRAGMENTS)







x  = AcGD3("AcGD3",[[18,1,1], [23,0,0]],  adduct="[M+H]+")
print x.printNist()


