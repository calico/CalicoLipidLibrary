from lipidRules import *

class AcGQ1b(SphingoLipid):

		neg_adduct_set = ["[M-2H]2-", "[M-3H]3-"]
		pos_adduct_set = ["[M+H]+", "[M+2H]2+"]

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = (MASS + ADDUCT[self.adduct])/abs(ADDUCT_CHARGE[self.adduct])

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 3, "precursor"])
					FRAGMENTS.append( [PREC-H2O, 73, "NL water"] )
					FRAGMENTS.append( [MW("C11H19NO9")+PROTON, 13, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O+PROTON, 403, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-2*H2O+PROTON, 639, "deoxy-N-Ac-Neuraminic Acid - water "] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9"), 47, "NL N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC - 2*MW("C11H19NO9")+H2O, 205, "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C6H12O6")+H2O, 4, "NL deoxy-N-Ac-Neuraminic Acid + Galactose + water "] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C6H12O6")+2*H2O, 16, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose "] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+2*H2O, 98, "NL deoxy-N-Ac-Neuraminic Acid + deoxyGalactose + N-Ac-Galactosamine "] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+3*H2O, 189, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H12O6")+4*H2O, 153, "NL 2x deoxy-N-Ac-Neuraminic Acid - 2x deoxyGalactose - N-Ac-Galactosamine"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+4*H2O, 81, "Ceramide - 2xwater"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O, 375, "Ceramide - water"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+6*H2O, 12, "Ceramide"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1])+H2O, 77, "sphingoid base "] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1]), 673, "sphingoid base - water"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1])+H2O-MW("CH2O"), 57, "sphingoid base - CH2O"] )
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 16, "sn2 fatty amide"])
					FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 588, "deoxy-N-Ac-glucosamine"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C11H19NO9") - 2*H2O + PROTON, 32, "deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3*H2O + PROTON, 16, "deoxy-N-Ac-galcatosamine + deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2*H2O + PROTON, 1000, "deoxy-N-Ac-glucosamine + deoxy-galactose"])

					FRAGMENTS.append( [PREC-MW("C6H12O6")-MW("C8H15NO6")+H2O, 77, "NL deoxygalactose + N-Ac-galactosamine"] )
					FRAGMENTS.append( [PREC-MW("C6H12O6"), 5, "NL galactose"] )

				if adduct == "[M+2H]2+":
					SINGLE = MASS + PROTON
					FRAGMENTS.append([PREC, 1, "precursor"])

					FRAGMENTS.append( [SINGLE - 2*MW("C11H19NO9")+H2O, 3, "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE- 2*MW("C11H19NO9")-MW("C6H12O6")+2*H2O, 1, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose "] )
					FRAGMENTS.append( [SINGLE- 2*MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+3*H2O, 11, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine"] )
					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C6H12O6")+3*H2O, 7, "NL N-Ac-Neuraminica Acid + 2x deoxy-N-Ac-Neuraminic Acid + Galactose "] )
					FRAGMENTS.append( [SINGLE - 4*MW("C11H19NO9")+3*H2O, 13, "NL N-Ac-Neuraminic Acid + 3x deoxy-N-acetyl-Neuraminic Acid"] )
					FRAGMENTS.append( [SINGLE- 3*MW("C11H19NO9")-MW("C6H12O6")-MW("C8H15NO6")+4*H2O, 35, "NL N-Ac-Neuraminica Acid + 2x deoxy-N-Ac-Neuraminic Acid + Galactose +N-Ac-Glucosamine"] )
					FRAGMENTS.append( [SINGLE- 4*MW("C11H19NO9")-MW("C6H12O6")-MW("C8H15NO6")+5*H2O, 51, "NL N-Ac-Neuraminica Acid + 3x deoxy-N-Ac-Neuraminic Acid + Galactose +N-Ac-Glucosamine"] )
					FRAGMENTS.append( [SINGLE- 4*MW("C11H19NO9")-2*MW("C6H12O6")-MW("C8H15NO6")+6*H2O, 58, "NL N-Ac-Neuraminica Acid + 3x deoxy-N-Ac-Neuraminic Acid + 2x Galactose +N-Ac-Glucosamine"] )
					FRAGMENTS.append( [SINGLE- 4*MW("C11H19NO9")-3*MW("C6H12O6")-MW("C8H15NO6")+7*H2O, 219, "Ceramide - water"] )
					FRAGMENTS.append( [SINGLE- 4*MW("C11H19NO9")-3*MW("C6H12O6")-MW("C8H15NO6")+6*H2O, 43, "Ceramide - 2x water"] )
					FRAGMENTS.append( [SINGLE- 4*MW("C11H19NO9")-3*MW("C6H12O6")-MW("C8H15NO6")+8*H2O - NL(self.chains[1]), 51, "sphingoid base "] )
					FRAGMENTS.append( [SINGLE- 4*MW("C11H19NO9")-3*MW("C6H12O6")-MW("C8H15NO6")+7*H2O - NL(self.chains[1]), 536, "sphingoid base - water"] )
					FRAGMENTS.append( [SINGLE- 4*MW("C11H19NO9")-3*MW("C6H12O6")-MW("C8H15NO6")+8*H2O - NL(self.chains[1]) - MW("CH2O"), 40, "sphingoid base - CH2O"] )

					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 5, "fatty amide"])
					FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 193, "deoxy-N-Ac-glucosamine"])
					FRAGMENTS.append([MW("C8H15NO6") - 2*H2O + PROTON, 109, "deoxy-N-Ac-glucosamine - water"])
					FRAGMENTS.append( [MW("C11H19NO9")-H2O+PROTON, 482, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-2*H2O+PROTON, 1000, "deoxy-N-Ac-Neuraminic Acid - water "] )

					FRAGMENTS.append([MW("C6H12O6") + MW("C11H19NO9") - 2*H2O + PROTON, 32, "deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2*H2O + PROTON, 92, "deoxy-N-Ac-glucosamine + deoxy-galactose"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3*H2O + PROTON, 41, "deoxy-N-Ac-galcatosamine + deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])

					FRAGMENTS.append([MW("C8H15NO6") + MW("C11H19NO9") - 2*H2O + PROTON, 2, "deoxy-N-Ac-galactosamine + deoxy-N-Ac-Neuraminic Acid"])

					FRAGMENTS.append([2*MW("C11H19NO9") - 2*H2O + PROTON, 17, "2x deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([2*MW("C11H19NO9") - 3*H2O + PROTON, 26, "2x deoxy-N-Ac-Neuraminic Acid - water"])




			elif adduct in NEG_ADDUCTS:  

				if adduct == "[M-2H]2-":
					SINGLE = MASS - PROTON
					FRAGMENTS.append([PREC, 21, "precursor"])
					FRAGMENTS.append( [PREC - H2O/2, 4, "NL water 2-"] )
					FRAGMENTS.append( [PREC - MW("CO2")/2, 6, "NL CO2 2-"] )
					FRAGMENTS.append( [PREC-(MW("C11H19NO9")-H2O)/2, 112, "NL deoxy-N-Ac-Neuraminic Acid 2-"] )
					FRAGMENTS.append( [PREC-(MW("C11H19NO9")-H2O+MW("CO2"))/2, 40, "NL deoxy-N-Ac-Neuraminic Acid + CO2 2-"] )
					FRAGMENTS.append( [PREC-2*(MW("C11H19NO9")-H2O)/2, 112, "NL 2x deoxy-N-Ac-Neuraminic Acid 2-"] )


					FRAGMENTS.append( [NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON, 5, "N-acyl-ethanolamine - H2O - H"])
					FRAGMENTS.append( [MW("C11H19NO9")-PROTON, 12, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O-PROTON, 1000, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-2*H2O-PROTON, 565, "2x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-2*H2O-PROTON, 132, "2x deoxy-N-Ac-Neuraminic Acid - CO2"] )
					FRAGMENTS.append( [MW("C11H19NO9")+MW("C11H19NO9")-MW("CO2")-3*H2O-PROTON, 32, "2x deoxy-N-Ac-Neuraminic Acid - CO2 - water"] )

					#FRAGMENTS.append( [SINGLE-MW("C11H19NO9")+H2O, 114, "NL deoxy-N-Ac-Neuraminic Acid"] )
					#FRAGMENTS.append( [SINGLE-MW("C11H19NO9"), 8, "NL N-Ac-Neuraminic Acid"] )
					#FRAGMENTS.append( [SINGLE-MW("C11H19NO9")-MW("CO2")+H2O, 24, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )
					FRAGMENTS.append( [SINGLE-2*MW("C11H19NO9")+2*H2O, 242, "NL 2x deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") -MW("CO2")+ 2 * H2O, 43, "NL 2x deoxy-N-Ac-Neuraminic Acid + CO2"])
					FRAGMENTS.append( [SINGLE-3*MW("C11H19NO9")+3*H2O, 634, "NL 3x deoxy-N-Ac-Neuraminic Acid"] )
					#FRAGMENTS.append( [SINGLE-3*MW("C11H19NO9")-MW("C6H10O5")+3*H2O, 4, "NL 3x deoxy-N-Ac-Neuraminic Acid + galactose"] )

					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 2, "Fatty Amide - NH3 - H"])
					FRAGMENTS.append([SINGLE - 4 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H10O5") + 5 * H2O, 18, "Ceramide"])
					FRAGMENTS.append([SINGLE - 4 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 4 * H2O, 4, "Glucosyl ceramide - water"])
					FRAGMENTS.append([SINGLE - 4 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 5* H2O, 10, "Glucosyl ceramide"])
					FRAGMENTS.append([SINGLE - 4 * MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H10O5") + 5 * H2O, 45, "Lactosyl ceramide"])

				if adduct == "[M-3H]3-":
					SINGLE = MASS - PROTON
					DOUBLE = (MASS - 2*PROTON)/2
					FRAGMENTS.append([PREC, 105, "precursor"])
					FRAGMENTS.append([SINGLE - 3 * MW("C11H19NO9") + 3 * H2O, 147, "NL 3x deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([SINGLE - 3 * MW("C11H19NO9") - MW("C6H12O6") + 3 * H2O, 147, "NL 3x deoxy-N-Ac-Neuraminic Acid"])

					FRAGMENTS.append([DOUBLE - 2 * (MW("C11H19NO9") - H2O) / 2, 297, "NL 2x deoxy-N-Ac-Neuraminic Acid 2-"])
					FRAGMENTS.append( [MW("C11H19NO9")-H2O-PROTON, 1000, "deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append([SINGLE - 4 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 4 * H2O, 55, "Glucosyl ceramide - water"])
					FRAGMENTS.append([MW("C6H12O6") - PROTON, 23, "Hexose"])
					FRAGMENTS.append([SINGLE, 23, "Single"])


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







# x  = AcGQ1b("AcGQ1b",[[18,1,1], [18,0,0]],  adduct="[M+2H]2+")
# print x.printNist()
