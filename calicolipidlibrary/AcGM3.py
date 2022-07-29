from lipidRules import *
#created in analagy to GcGM3 frag patterns from DO plasma 
class AcGM3(SphingoLipid):


		pos_adduct_set = ["[M+H]+"]
		neg_adduct_set = ["[M-H]-"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
				#couldn't modify -- no spectrum
				if adduct == "[M+H]+":     
					FRAGMENTS.append([PREC, 10, "precursor"])
					FRAGMENTS.append([PREC-H2O, 25, "NL water"])
					FRAGMENTS.append( [PREC-MW("C11H19NO9"), 65, "NL N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C6H10O5"), 114, "NL N-Ac-Neuraminic Acid + galactose"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-2*MW("C6H10O5"), 409, "Ceramide"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-2*MW("C6H10O5")-H2O, 97, "Ceramide - water"] )
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-2*MW("C6H10O5")-NL(self.chains[1])+H2O, 116, "LCB"])
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-2*MW("C6H10O5")-NL(self.chains[1]), 1000, "LCB - water"])
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-2*MW("C6H10O5")-NL(self.chains[1])+H2O - MW("CH2O"), 90, "LCB - CH2O"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 22, "fatty amide"])


					FRAGMENTS.append( [MW("C11H19NO9")-H2O+PROTON, 38, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-2*H2O+PROTON, 101, "deoxy-N-Ac-Neuraminic Acid - water"] )



			elif adduct in NEG_ADDUCTS:
				if adduct == "[M-H]-":
					FRAGMENTS.append( [PREC, 1000, "Precursor"] )
					FRAGMENTS.append( [PREC - H2O, 4, "NL water"] )
					FRAGMENTS.append( [PREC - MW("CO2"), 10, "NL CO2"] )

					FRAGMENTS.append( [MW("C11H19NO9")-PROTON, 2, "N-Ac-Neuraminic Acid "] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O-PROTON, 837, "deoxy-N-Ac-Neuraminic Acid"] )

					FRAGMENTS.append( [MW("C6H12O6")-H2O-PROTON, 8, "deoxyglucose"] )
					FRAGMENTS.append( [MW("C6H12O6")-PROTON, 7, "glucose"] )

					FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 14, "NL - N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C6H10O5")+H2O, 4, "NL N-Ac-Neuraminic Acid - galactose"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-2*MW("C6H10O5")+H2O, 21, "Ceramide"] )


			
			return(FRAGMENTS)

# calicolipidlibrary.print_spectrum("AcGM3", [[18,1,1], [24,1,0]], "[M-H]-")

# x  = AcGM3("AcGM3",[[18,1,1], [24,1,0]],  adduct="[M-H]-")
# print x.printNist()
#
# x  = AcGM3("AcGM3",[[18,1,1], [23,0,0]],  adduct="[M+H]+")
# print x.printNist()

