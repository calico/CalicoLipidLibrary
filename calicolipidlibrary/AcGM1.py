from lipidRules import *
#created from DO plasma frag patterns
class AcGM1(SphingoLipid):

		pos_adduct_set = ["[M+H]+"]
		neg_adduct_set = ["[M-H]-"]

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 1, "precursor"])

					FRAGMENTS.append( [PREC-H2O, 26, "NL water"] )
					FRAGMENTS.append( [PREC - MW("C6H12O6"), 3, "NL Galactose"] )
					FRAGMENTS.append( [PREC - MW("C6H12O6")-MW("C8H15NO6")+H2O, 64, "NL Galactose + GalNAc"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9"), 40, "NL N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC - MW("C6H12O6")-MW("C8H15NO6")-MW("C11H19NO9")+2*H2O, 82, "NL Galactose + GalNAc + N-Ac_Neuraminic Acid"] )

					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")+H2O, 136, "NL N-Ac-Neuraminic Acid - N-Ac-Glucosamine - glucose - galactose"] )
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H10O5")-NL(self.chains[1])+H2O*2, 91, "Sphingoid base"])
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H10O5")-NL(self.chains[1])+H2O, 807, "Sphingoid base - water"])
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H10O5")-NL(self.chains[1])+H2O*2-MW("CH2O"), 72, "Sphingoid base - CH2O"])
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H10O5")+H2O*2, 6, "Ceramide"])
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H10O5")+H2O, 397, "Ceramide - H2O"])
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H10O5"), 81, "Ceramide - 2x H2O"])

					FRAGMENTS.append([MW("C8H15NO6") + PROTON, 8, "N-Ac-glucosamine"])
					FRAGMENTS.append( [MW("C8H15NO6")-H2O+PROTON, 1000, "N-Ac-glucosamine - water"] )
					FRAGMENTS.append( [MW("C8H15NO6")-2*H2O+PROTON, 853, "N-Ac-glucosamine - 2*water"] )
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 91, "fatty amide"])
					FRAGMENTS.append( [MW("C11H19NO9")-H2O+PROTON, 64, "deoxy-N-Ac-Neuraminic acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C6H12O6")-2*H2O+PROTON, 4, "deoxy-N-Ac-Neuraminic acid + deoxygalactose"] )
					FRAGMENTS.append( [MW("C8H15NO6")+MW("C6H12O6")-2*H2O+PROTON, 980, "deoxy-galNAc + deoxygalactose"] )


			
			elif adduct in NEG_ADDUCTS:  
				if adduct == "[M+Cl]-" or adduct == "[M+FA]-":
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON, 1000, "NL adduct"] )
				if adduct == "[M-H]-":
					FRAGMENTS.append( [PREC, 1000, "Precursor"] )
					FRAGMENTS.append( [PREC - H2O, 5, "NL water"] )
					FRAGMENTS.append( [PREC - MW("CO2"), 4, "NL CO2"] )
					FRAGMENTS.append( [MW("C11H19NO9")-PROTON, 10, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O-PROTON, 296, "deoxy-N-Ac-Neuraminic Acid"] )

					FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 20, "NL deoxy N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C6H10O5")+H2O, 1, "NL deoxy N-Ac-Neuraminic Acid + deoxyGalactose"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C6H10O5"), 2, "NL deoxy N-Ac-Neuraminic Acid + Galactose"] )

					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H10O5")+2*H2O, 43, "NL N-Ac-Neuraminic Acid + Galactose + N Ac Glucosamine"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")+2*H2O, 16, "Glucosyl Ceramide"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")+1*H2O, 7, "Glucosyl Ceramide - H2O"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H10O5")+2*H2O, 25, "Ceramide"] )


			return(FRAGMENTS)







# x  = AcGM1("AcGM1",[[18,1,1], [18,0,0]],  adduct="[M-H]-")
# print x.printNist()
#
#
# x  = AcGM1("AcGM1",[[18,1,1], [18,0,0]],  adduct="[M+H]+")
# print x.printNist()