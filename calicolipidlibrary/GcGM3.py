from lipidRules import *
#created from DO plasma frag patterns
class GcGM3(SphingoLipid):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append( [MW("C11H19NO10")-H2O+PROTON, 1000, "N-Ac-Neuraminic Acid - water + H+"] )				  #??
					FRAGMENTS.append( [PREC-H2O, 1000, "pre - water"] )					
				
					FRAGMENTS.append( [PREC-MW("C11H19NO10"), 1000, "pre - N-Gc-Neuraminic Acid"] )					
					FRAGMENTS.append( [PREC-MW("C11H19NO10")-MW("C6H10O5"), 1000, "pre - N-Gc-Neuraminic Acid - glucose"] )									
					FRAGMENTS.append( [PREC-MW("C11H19NO10")-2*MW("C6H10O5"), 1000, "pre - N-Gc-Neuraminic Acid - glucose - galactose"] )													
					FRAGMENTS.append( [PREC-MW("C11H19NO10")-2*MW("C6H10O5")-H2O, 1000, "pre - N-Gc-Neuraminic Acid - glucose - galactose - water"] )													
					FRAGMENTS.append( [PREC - MW("C11H19NO10")-2*MW("C6H10O5")-NL(self.chains[1])+H2O, 1000, "pre- (head + sn2)"])

					FRAGMENTS.append( [PREC - MW("C11H19NO10")-2*MW("C6H10O5")-NL(self.chains[1]), 1000, "pre - (head + sn2 -water )"])
					FRAGMENTS.append( [PREC - MW("C11H19NO10")-2*MW("C6H10O5")-NL(self.chains[1])-MW("C"), 1000, "pre - (head + sn2 -CH2O )"])					
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 1000, "sn2 fatty amide"])
				
			 

				if adduct == "[M+Na]+" or adduct == "[M+K]+":
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append( [PREC-MW("C11H19NO10")+H2O, 1000, "pre - N-Gc-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO10")-MW("C6H10O5")+H2O, 1000, "pre - N-Gc-Neuraminic Acid - glucose"] )									
					FRAGMENTS.append( [PREC-MW("C11H19NO10")-2*MW("C6H10O5")+H2O, 1000, "pre - N-Gc-Neuraminic Acid - glucose - galactose"] )													
					FRAGMENTS.append( [PREC-MW("C11H19NO10")-2*MW("C6H10O5"), 1000, "pre - N-Gc-Neuraminic Acid - glucose - galactose - water"] )													
					#FRAGMENTS.append( [PREC - MW("C11H19NO10")-2*MW("C6H10O5")-NL(self.chains[1])+H2O, 1000, "pre- (head + sn2)"])

					#FRAGMENTS.append( [NL(self.chains[1])-H2O+PROTON, 1000, "sn2 acylium ion"])
					#FRAGMENTS.append( [NL(self.chains[1])-H2O-H2O+PROTON, 1000, "sn2 acylium ion - H2O"])
	
			
			elif adduct in NEG_ADDUCTS:
				if adduct == "[M-H]-": #by analogy from AcGM3 standard
					FRAGMENTS.append([PREC, 1000, "Precursor"])
					FRAGMENTS.append([PREC - H2O, 4, "NL water"])
					FRAGMENTS.append([PREC - MW("CO2"), 10, "NL CO2"])

					FRAGMENTS.append([MW("C11H19NO10") - PROTON, 2, "N-Ac-Neuraminic Acid "])
					FRAGMENTS.append([MW("C11H19NO10") - H2O - PROTON, 837, "deoxy-N-Ac-Neuraminic Acid"])

					FRAGMENTS.append([MW("C6H12O6") - H2O - PROTON, 8, "deoxyglucose"])
					FRAGMENTS.append([MW("C6H12O6") - PROTON, 7, "glucose"])

					FRAGMENTS.append([PREC - MW("C11H19NO10") + H2O, 14, "NL - N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([PREC - MW("C11H19NO10") - MW("C6H10O5") + H2O, 4, "NL N-Ac-Neuraminic Acid - galactose"])
					FRAGMENTS.append([PREC - MW("C11H19NO10") - 2 * MW("C6H10O5") + H2O, 21, "Ceramide"])


			return(FRAGMENTS)    




# x  = GcGM3("GcGM3",[[18,1,1], [24,1,0]],  adduct="[M-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("GcGM3", [[18,1,1], [24,1,0]], "[M-H]-")