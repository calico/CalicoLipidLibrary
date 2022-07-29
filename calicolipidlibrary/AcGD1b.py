from lipidRules import *

class AcGD1b(SphingoLipid):

		neg_adduct_set = ["[M-H]-", "[M-2H]2-"]
		pos_adduct_set = Lipid.pos_adduct_set + ["[M+H]+", "[M+2H]2+", "[M+Na]+"]

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = (MASS + ADDUCT[self.adduct])/abs(ADDUCT_CHARGE[self.adduct])

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 1, "precursor"])
					FRAGMENTS.append( [PREC-H2O, 30, "NL water"] )
					#FRAGMENTS.append( [PREC-MW("C6H12O6"), 5, "NL galactose"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9"), 23, "NL N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C6H12O6")-MW("C8H15NO6")+H2O, 41, "NL deoxygalactose + N-Ac-galactosamine"] )

					FRAGMENTS.append( [PREC - 2*MW("C11H19NO9")+H2O, 125, "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid"] )
					#FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C6H12O6")+H2O, 4, "NL deoxy-N-Ac-Neuraminic Acid + Galactose + water "] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+2*H2O, 65, "NL deoxy-N-Ac-Neuraminic Acid + deoxyGalactose + N-Ac-Galactosamine "] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C6H12O6")+2*H2O, 12, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose "] )

					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+3*H2O, 151, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H12O6")+4*H2O, 126, "NL 2x deoxy-N-Ac-Neuraminic Acid + 2x deoxyGalactose + N-Ac-Galactosamine"] )

					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+4*H2O, 71, "Ceramide - 2xwater"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O, 323, "Ceramide - water"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+6*H2O, 4, "Ceramide"] )

					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1])+H2O, 78, "sphingoid base "] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1]), 726, "sphingoid base - water"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1])+H2O-MW("CH2O"), 63, "sphingoid base - CH2O"] )
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 17, "sn2 fatty amide"])

					FRAGMENTS.append( [MW("C11H19NO9")+PROTON, 12, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O+PROTON, 413, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-2*H2O+PROTON, 667, "deoxy-N-Ac-Neuraminic Acid - water "] )

					FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 851, "deoxy-N-Ac-glucosamine"])
					FRAGMENTS.append([MW("C8H15NO6") - 2*H2O + PROTON, 660, "deoxy-N-Ac-glucosamine - water"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C11H19NO9") - 2*H2O + PROTON, 26, "deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3*H2O + PROTON, 3, "deoxy-N-Ac-galcatosamine + deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2*H2O + PROTON, 1000, "deoxy-N-Ac-glucosamine + deoxy-galactose"])

			if adduct == "[M+2H]2+":
				SINGLE = MASS + PROTON

				FRAGMENTS.append([PREC, 1, "precursor"])
				#FRAGMENTS.append([PREC - H2O, 30, "NL water"])
				#FRAGMENTS.append( [PREC-MW("C6H12O6")/2, 5, "NL galactose"] )
				FRAGMENTS.append([SINGLE - MW("C11H19NO9"), 10, "NL N-Ac-Neuraminic Acid"])
				FRAGMENTS.append([SINGLE - MW("C6H12O6") - MW("C8H15NO6") + H2O, 8, "NL deoxygalactose + N-Ac-galactosamine"])

				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") + H2O, 41, "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid"])
				# FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C6H12O6")+H2O, 4, "NL deoxy-N-Ac-Neuraminic Acid + Galactose + water "] )
				FRAGMENTS.append([SINGLE - MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 2 * H2O, 62, "NL deoxy-N-Ac-Neuraminic Acid + deoxyGalactose + N-Ac-Galactosamine "])
				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C6H12O6") + 2 * H2O, 64, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose "])

				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 3 * H2O, 108, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine"])
				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H12O6") + 4 * H2O, 115, "NL 2x deoxy-N-Ac-Neuraminic Acid + 2x deoxyGalactose + N-Ac-Galactosamine"])

				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 4 * H2O, 88, "Ceramide - 2xwater"])
				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O, 341, "Ceramide - water"])
				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 6 * H2O, 2, "Ceramide"])

				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]) + H2O, 105, "sphingoid base "])
				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]), 1000, "sphingoid base - water"])
				FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]) + H2O - MW("CH2O"), 80, "sphingoid base - CH2O"])
				FRAGMENTS.append([NL(self.chains[1]) - H2O + MW("NH3") + PROTON, 14, "sn2 fatty amide"])

				#FRAGMENTS.append([MW("C11H19NO9") + PROTON, 12, "N-Ac-Neuraminic Acid"])
				FRAGMENTS.append([MW("C11H19NO9") - H2O + PROTON, 284, "deoxy-N-Ac-Neuraminic Acid"])
				FRAGMENTS.append([MW("C11H19NO9") - 2 * H2O + PROTON, 898, "deoxy-N-Ac-Neuraminic Acid - water "])

				FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 635, "deoxy-N-Ac-glucosamine"])
				FRAGMENTS.append([MW("C8H15NO6") - 2 * H2O + PROTON, 461, "deoxy-N-Ac-glucosamine - water"])
				FRAGMENTS.append([MW("C6H12O6") + MW("C11H19NO9") - 2 * H2O + PROTON, 3, "deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
				#FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3 * H2O + PROTON, 10, "deoxy-N-Ac-galcatosamine + deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
				FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + PROTON, 417, "deoxy-N-Ac-glucosamine + deoxy-galactose"])

			if adduct == "[M+Na]+" or adduct == "[M+K]+" or adduct == "[M+Li]+":
				FRAGMENTS.append([PREC, 98, "precursor"])
				FRAGMENTS.append([PREC - H2O, 5, "NL water"])
				FRAGMENTS.append([PREC - MW("CO2"), 4, "NL CO2"])

				FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 65, "NL deoxy-N-Ac-Neuraminic Acid"] )
				FRAGMENTS.append( [PREC-MW("C11H19NO9"), 2, "NL N-Ac-Neuraminic Acid"] )
				FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") + 2*H2O, 1000, "NL 2x deoxy-N-Ac-Neuraminic Acid"])
				FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C6H12O6") + 3*H2O, 20, "NL 2x deoxy-N-Ac-Neuraminic Acid + deoxygalactose"])

				FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 4 * H2O, 339,
								  "NL 2x deoxy-N-Ac-Neuraminic Acid + deoxy-Galactose + deoxy-N-Ac-Galactosamine"])
				FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H12O6") + 5 * H2O, 43,
								  "NL 2x deoxy-N-Ac-Neuraminic Acid + 2x deoxyGalactose + deoxy-N-Ac-Galactosamine"])
				FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - ADDUCT[adduct] + PROTON, 1,
								  "Ceramide - water - Adduct"])
				FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]) - ADDUCT[adduct] + PROTON, 5,
								  "sphingoid base - water - ADDUCT"])
				FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 10, "deoxy-N-Ac-glucosamine"])
				FRAGMENTS.append([MW("C8H15NO6") - 2 * H2O + PROTON, 8, "deoxy-N-Ac-glucosamine - water"])
				FRAGMENTS.append([MW("C11H19NO9") - H2O + PROTON, 4, "deoxy-N-Ac-Neuraminic Acid"])
				FRAGMENTS.append([MW("C11H19NO9") - 2 * H2O + PROTON, 7, "deoxy-N-Ac-Neuraminic Acid - water "])
				FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + PROTON, 12, "deoxy-N-Ac-glucosamine + deoxygalactose"])
				FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + ADDUCT[adduct], 81, "deoxy-N-Ac-glucosamine + deoxygalactose + ADDUCT"])


			elif adduct in NEG_ADDUCTS:  
				if adduct == "[M+Cl]-" or adduct == "[M+FA]-":
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON, 1000, "NL adduct"] )
				if adduct == "[M-H]-":
					FRAGMENTS.append([PREC, 304, "precursor"])
					FRAGMENTS.append( [PREC - H2O, 10, "NL water"] )
					FRAGMENTS.append( [PREC - MW("CO2"), 49, "NL CO2"] )
					FRAGMENTS.append( [MW("C11H19NO9")-PROTON, 4, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O-PROTON, 432, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-2*H2O-PROTON, 164, "2x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-2*H2O-PROTON, 54, "2x deoxy-N-Ac-Neuraminic Acid - CO2"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-3*H2O-PROTON, 13, "2x deoxy-N-Ac-Neuraminic Acid - CO2 - water"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 1000, "NL deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9"), 4, "NL N-Ac-Neuraminic Acid"] )

					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("CO2")+H2O, 2, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )

					FRAGMENTS.append( [PREC-2*MW("C11H19NO9")+2*H2O, 15, "NL 2x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON, 8, "N-acyl-ethanolamine - H2O - H"])
					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 2, "Fatty Amide - NH3 - H"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H10O5") + 3 * H2O, 24, "Ceramide"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 2 * H2O, 7, "Glucosyl ceramide - water"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 3* H2O, 17, "Glucosyl ceramide"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") -  MW("C6H10O5") + 3 * H2O, 48, "Lactosyl ceramide"])

				if adduct == "[M-2H]2-":
					SINGLE = MASS - PROTON
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append( [PREC - H2O/2, 11, "NL water 2-"] )
					FRAGMENTS.append( [PREC - MW("CO2")/2, 55, "NL CO2 2-"] )
					FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON, 16, "N-acyl-ethanolamine - H2O - H"])
					FRAGMENTS.append( [MW("C11H19NO9")-PROTON, 19, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O-PROTON, 612, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-2*H2O-PROTON, 326, "2x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-2*H2O-PROTON, 11, "2x deoxy-N-Ac-Neuraminic Acid - CO2"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-3*H2O-PROTON, 4, "2x deoxy-N-Ac-Neuraminic Acid - CO2 - water"] )

					FRAGMENTS.append( [SINGLE-MW("C11H19NO9")+H2O, 153, "NL deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE-MW("C11H19NO9"), 3, "NL N-Ac-Neuraminic Acid"] )
					#FRAGMENTS.append( [SINGLE-MW("C11H19NO9")-MW("CO2")+H2O, 5, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )
					FRAGMENTS.append( [SINGLE-2*MW("C11H19NO9")+2*H2O, 132, "NL 2x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE-2*MW("C11H19NO9")+2*H2O-MW("C6H10O5"), 3, "NL 2x deoxy-N-Ac-Neuraminic Acid + deoxygalactose"] )


					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 2, "Fatty Amide - NH3 - H"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H10O5") + 3 * H2O, 49, "Ceramide"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 2 * H2O, 15, "Glucosyl ceramide - water"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 3* H2O, 30, "Glucosyl ceramide"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") -  MW("C6H10O5") + 3 * H2O, 139, "Lactosyl ceramide"])


			return(FRAGMENTS)






x  = AcGD1b("AcGD1b",[[18,1,1], [18,0,0]],  adduct="[M+H]+")
print x.printNist()


x  = AcGD1b("AcGD1b",[[18,1,1], [18,0,0]],  adduct="[M+Na]+")
print x.printNist()
