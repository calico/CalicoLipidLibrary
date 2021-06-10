from lipidRules import *

class AcGT1b(SphingoLipid):

		neg_adduct_set = ["[M-2H]2-"]
		pos_adduct_set = ["[M+H]+", "[M+2H]2+", "[M+Na]+"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = (MASS + ADDUCT[self.adduct])/abs(ADDUCT_CHARGE[self.adduct])

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 3, "precursor"])
					FRAGMENTS.append([PREC - H2O, 73, "NL water"])

					FRAGMENTS.append([MW("C11H19NO9") - H2O + PROTON, 514, "deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C11H19NO9") - 2 * H2O + PROTON, 18, "deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + PROTON, 112, "deoxy-N-Ac-glucosamine + deoxy-galactose"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3 * H2O + PROTON, 204, "deoxy-N-Ac-galcatosamine + deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])

					FRAGMENTS.append([MW("C11H19NO9") - 2 * H2O + PROTON, 1000, "deoxy-N-Ac-Neuraminic Acid - water "])
					FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 425, "deoxy-N-Ac-glucosamine"])

					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") + H2O, 3, "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid"])
					FRAGMENTS.append([PREC - MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 2 * H2O, 6, "NL deoxy-N-Ac-Neuraminic Acid + deoxyGalactose + N-Ac-Galactosamine "])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 3 * H2O, 6, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine"])

					FRAGMENTS.append([PREC - 3 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O, 3, "Ceramide "])
					FRAGMENTS.append([PREC - 3 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 6 * H2O, 79, "Ceramide - water"])
					FRAGMENTS.append([PREC - 3 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 6 * H2O - NL(self.chains[1]) + H2O, 3, "sphingoid base"])
					FRAGMENTS.append([PREC - 3 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]) + H2O, 270, "sphingoid base - water "])



				if adduct == "[M+2H]2+":
					SINGLE = MASS + PROTON
					FRAGMENTS.append([PREC, 1, "precursor"])

					FRAGMENTS.append( [SINGLE-MW("C11H19NO9"), 2, "NL N-Ac-Neuraminic Acid"] )

					FRAGMENTS.append( [SINGLE - 2*MW("C11H19NO9")+H2O, 9, "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE-MW("C11H19NO9")-MW("C6H12O6")+H2O, 6, "NL deoxy-N-Ac-Neuraminic Acid + deoxygalactose"] )
					FRAGMENTS.append( [SINGLE- 2*MW("C11H19NO9")-MW("C6H12O6")+2*H2O, 17, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose "] )
					FRAGMENTS.append( [SINGLE - 3*MW("C11H19NO9")+2*H2O, 14, "NL N-Ac-Neuraminic Acid + 2x deoxy-N-acetyl-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE-MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+2*H2O, 12, "NL deoxy-N-Ac-Neuraminic Acid + deoxyGalactose + N-Ac-Galactosamine "] )
					FRAGMENTS.append( [SINGLE- 2*MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+3*H2O, 36, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine"] )
					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+4*H2O, 59, "NL N-Ac-Neuaraminic Acid + 3x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine"] )
					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H12O6")+5*H2O, 68, "NL N-Ac-Neuaraminic Acid + 3x deoxy-N-Ac-Neuraminic Acid + 2*Galactose + N-Ac-Galactosamine"] )

					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O, 54, "Ceramide - 2xwater"] )
					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+6*H2O, 223, "Ceramide - water"] )
					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+7*H2O, 2, "Ceramide"] )

					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+6*H2O - NL(self.chains[1])+H2O, 64, "sphingoid base "] )
					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1])+H2O, 613, "sphingoid base - water "] )
					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+6*H2O - NL(self.chains[1])+H2O-MW("CH2O"), 49, "sphingoid base - CH2O"] )
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 11, "sn2 fatty amide"])

					FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 286, "deoxy-N-Ac-glucosamine"])
					FRAGMENTS.append([MW("C8H15NO6") - 2*H2O + PROTON, 166, "deoxy-N-Ac-glucosamine - water"])

					FRAGMENTS.append( [MW("C11H19NO9")+PROTON, 5, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O+PROTON, 404, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-2*H2O+PROTON, 1000, "deoxy-N-Ac-Neuraminic Acid - water "] )


					FRAGMENTS.append([MW("C6H12O6") + MW("C11H19NO9") - 2*H2O + PROTON, 22, "deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3*H2O + PROTON, 69, "deoxy-N-Ac-galcatosamine + deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2*H2O + PROTON, 86, "deoxy-N-Ac-glucosamine + deoxy-galactose"])


			if adduct == "[M+Na]+" or adduct == "[M+K]+":
				FRAGMENTS.append([PREC, 56, "precursor"])

				FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 135, "NL deoxy-N-Ac-Neuraminic Acid"] )
				FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") + 2*H2O, 319, "NL 2x deoxy-N-Ac-Neuraminic Acid"])
				FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") + 2*H2O - MW("CO2"), 38, "NL 2x deoxy-N-Ac-Neuraminic Acid + CO2"])
				FRAGMENTS.append([PREC - 3 * MW("C11H19NO9") + 3*H2O, 1000, "NL 3x deoxy-N-Ac-Neuraminic Acid"])

				FRAGMENTS.append([PREC - 3 * MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 5 * H2O, 443,
								  "NL 3x deoxy-N-Ac-Neuraminic Acid + deoxy-Galactose + deoxy-N-Ac-Galactosamine"])
				FRAGMENTS.append([PREC - 3 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H12O6") + 6 * H2O, 37,
								  "NL 2x deoxy-N-Ac-Neuraminic Acid + 2x deoxyGalactose + deoxy-N-Ac-Galactosamine"])

				FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + ADDUCT[adduct], 91, "deoxy-N-Ac-glucosamine + deoxygalactose + ADDUCT"])
				FRAGMENTS.append([MW("C11H19NO9") - H2O + PROTON, 1, "deoxy-N-Ac-Neuraminic Acid"])
				FRAGMENTS.append([MW("C11H19NO9") - 2 * H2O + PROTON, 3, "deoxy-N-Ac-Neuraminic Acid - water "])


			elif adduct in NEG_ADDUCTS:  

				if adduct == "[M-2H]2-":
					SINGLE = MASS - PROTON
					FRAGMENTS.append([PREC, 65, "precursor"])
					FRAGMENTS.append( [PREC - H2O/2, 24, "NL water 2-"] )
					FRAGMENTS.append( [PREC - MW("CO2")/2, 45, "NL CO2 2-"] )
					FRAGMENTS.append( [PREC-(MW("C11H19NO9")-H2O)/2, 670, "NL deoxy-N-Ac-Neuraminic Acid 2-"] )

					FRAGMENTS.append( [NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON, 9, "N-acyl-ethanolamine - H2O - H"])
					FRAGMENTS.append( [MW("C11H19NO9")-PROTON, 6, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O-PROTON, 1000, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-2*H2O-PROTON, 107, "2x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-2*H2O-PROTON, 32, "2x deoxy-N-Ac-Neuraminic Acid - CO2"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-3*H2O-PROTON, 7, "2x deoxy-N-Ac-Neuraminic Acid - CO2 - water"] )

					FRAGMENTS.append( [SINGLE-MW("C11H19NO9")+H2O, 114, "NL deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE-MW("C11H19NO9"), 8, "NL N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE-MW("C11H19NO9")-MW("CO2")+H2O, 24, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )
					FRAGMENTS.append( [SINGLE-2*MW("C11H19NO9")+2*H2O, 450, "NL 2x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE-3*MW("C11H19NO9")+3*H2O, 20, "NL 3x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE-3*MW("C11H19NO9")-MW("C6H10O5")+3*H2O, 4, "NL 3x deoxy-N-Ac-Neuraminic Acid + galactose"] )

					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 2, "Fatty Amide - NH3 - H"])
					FRAGMENTS.append([SINGLE - 3 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H10O5") + 4 * H2O, 21, "Ceramide"])
					FRAGMENTS.append([SINGLE - 3 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 3 * H2O, 6, "Glucosyl ceramide - water"])
					FRAGMENTS.append([SINGLE - 3 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 4* H2O, 13, "Glucosyl ceramide"])
					FRAGMENTS.append([SINGLE - 3 * MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H10O5") + 4 * H2O, 50, "Lactosyl ceramide"])


			return(FRAGMENTS)


		@classmethod
		def generateLibrary(cls,target=None, mode="pos"):
			if target: handle = open(target, 'a+')
			if mode == "pos":  
				adduct_set = cls.pos_adduct_set
			elif mode =="neg":
				adduct_set = cls.neg_adduct_set
			parent = cls.__bases__[0]
			class_name = cls.__name__ 
			for c in parent.chain_sets:
				for adduct in adduct_set:
					x  = cls(class_name, c, adduct=adduct)						
					content = x.printNist()
					if target: handle.write(content)
					else: sys.stdout.write(content)
			if target: handle.close()




x  = AcGT1b("AcGT1b",[[18,1,1], [18,0,0]],  adduct="[M+H]+")
print x.printNist()



x  = AcGT1b("AcGT1b",[[18,1,1], [18,0,0]],  adduct="[M+2H]2+")
print x.printNist()



x  = AcGT1b("AcGT1b",[[18,1,1], [18,0,0]],  adduct="[M+Na]+")
print x.printNist()

