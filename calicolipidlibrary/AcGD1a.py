from lipidRules import *

class AcGD1a(SphingoLipid):

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
					FRAGMENTS.append( [PREC-H2O, 3, "NL water"] )
					
					FRAGMENTS.append( [PREC-MW("C11H19NO9"), 12, "NL N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC - 2*MW("C11H19NO9")+H2O, 10, "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C6H12O6")+H2O, 2, "NL deoxy-N-Ac-Neuraminic Acid + Galactose + water "] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C6H12O6")+2*H2O, 2, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose "] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+2*H2O, 27, "NL deoxy-N-Ac-Neuraminic Acid + deoxyGalactose + N-Ac-Galactosamine "] )

					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H12O6")+3*H2O, 24, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H12O6")+4*H2O, 32, "NL 2x deoxy-N-Ac-Neuraminic Acid + 2x deoxyGalactose + N-Ac-Galactosamine"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+4*H2O, 18, "Ceramide - 2xwater"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O, 94, "Ceramide - water"] )
					#FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+6*H2O, 3, "Ceramide"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1])+H2O, 21, "sphingoid base "] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1]), 203, "sphingoid base - water"] )
					FRAGMENTS.append( [PREC- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+5*H2O - NL(self.chains[1])+H2O-MW("CH2O"), 17, "sphingoid base - CH2O"] )

					FRAGMENTS.append([MW("C11H19NO9") + PROTON, 10, "N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9") - H2O + PROTON, 534, "deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9") - 2 * H2O + PROTON, 1000, "deoxy-N-Ac-Neuraminic Acid - water "])
					FRAGMENTS.append([MW("C11H19NO9") - 3 * H2O + PROTON, 74, "deoxy-N-Ac-Neuraminic Acid - 2 water "])

					#FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 6, "sn2 fatty amide"])
					FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 527, "N-Ac-glucosamine - water"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C11H19NO9") - 2*H2O + PROTON, 38, "deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3*H2O + PROTON, 243, "deoxy-N-Ac-galcatosamine + deoxygalactose + dexoy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2*H2O + PROTON, 130, "deoxy-N-Ac-glucosamine + deoxy-galactose"])

				if adduct == "[M+2H]2+":
					SINGLE = MASS + PROTON

					FRAGMENTS.append([PREC, 2, "precursor"])

					#FRAGMENTS.append([SINGLE - MW("C11H19NO9"), 12, "NL N-Ac-Neuraminic Acid"])
					#FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") + H2O, 1, "NL N-Ac-Neuraminic Acid + deoxy-N-acetyl-Neuraminic Acid"])
					FRAGMENTS.append([SINGLE - MW("C11H19NO9") - MW("C6H12O6") + H2O, 14, "NL deoxy-N-Ac-Neuraminic Acid + Galactose + water "])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C6H12O6") + 2 * H2O, 33, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose "])
					FRAGMENTS.append([SINGLE - MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 2 * H2O, 34, "NL deoxy-N-Ac-Neuraminic Acid + deoxyGalactose + N-Ac-Galactosamine "])

					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 3 * H2O, 53, "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose + N-Ac-Galactosamine"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H12O6") + 4 * H2O, 100, "NL 2x deoxy-N-Ac-Neuraminic Acid + 2x deoxyGalactose + N-Ac-Galactosamine"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 4 * H2O, 50, "Ceramide - 2xwater"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O, 317, "Ceramide - water"])
					#FRAGMENTS.append( [SINGLE- 2*MW("C11H19NO9")-MW("C8H15NO6")-3*MW("C6H12O6")+6*H2O, 3, "Ceramide"] )
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]) + H2O, 57, "sphingoid base "])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]), 945, "sphingoid base - water"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]) + H2O - MW("CH2O"), 45, "sphingoid base - CH2O"])

					#FRAGMENTS.append([MW("C11H19NO9") + PROTON, 10, "N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9") - H2O + PROTON, 318, "deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9") - 2 * H2O + PROTON, 1000, "deoxy-N-Ac-Neuraminic Acid - water "])
					FRAGMENTS.append([MW("C11H19NO9") - 3 * H2O + PROTON, 39, "deoxy-N-Ac-Neuraminic Acid - 2 water "])

					#FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 6, "sn2 fatty amide"])
					FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 190, "N-Ac-glucosamine - water"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C11H19NO9") - 2 * H2O + PROTON, 2, "deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3 * H2O + PROTON, 36, "deoxy-N-Ac-galcatosamine + deoxygalactose + dexoy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + PROTON, 87, "deoxy-N-Ac-glucosamine + deoxy-galactose"])

				if adduct == "[M+Na]+" or adduct == "[M+K]+":
					FRAGMENTS.append([PREC, 166, "precursor"])
					FRAGMENTS.append([PREC - MW("CO2"), 4, "NL CO2"])

					FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 382, "NL deoxy-N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O-MW("CO2"), 90, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )

					FRAGMENTS.append( [PREC-MW("C11H19NO9"), 6, "NL N-Ac-Neuraminic Acid"] )

					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") + 2*H2O, 1000, "NL 2x deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([PREC - MW("C11H19NO9") - MW("C6H12O6") + H2O - ADDUCT[adduct] + PROTON, 4, "NL N-Ac-Neuraminic Acid + Galactose + ADDUCT "])
					FRAGMENTS.append([PREC - MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 3 * H2O -MW("CO2"), 28,
									  "NL deoxy-N-Ac-Neuraminic Acid + deoxyGalactose + N-Ac-Galactosamine -CO2"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C6H12O6") + 3 * H2O, 44,
									  "NL 2x deoxy-N-Ac-Neuraminic Acid + Galactose"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H12O6") + 4 * H2O, 435,
									  "NL 2x deoxy-N-Ac-Neuraminic Acid + deoxy-Galactose + deoxy-N-Ac-Galactosamine"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H12O6") + 5 * H2O, 51,
									  "NL 2x deoxy-N-Ac-Neuraminic Acid + 2x deoxyGalactose + deoxy-N-Ac-Galactosamine"])

					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 6 * H2O, 3,
									  "Ceramide + adduct"])


					#FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]) + H2O - ADDUCT[adduct] + PROTON, 12,
					#				  "sphingoid base - adduct "])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]) - ADDUCT[adduct] + PROTON, 2,
									  "sphingoid base - water - ADDUCT"])
					#FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H12O6") + 5 * H2O - NL(self.chains[1]) + H2O - MW("CH2O") - ADDUCT[adduct] + PROTON, 9,
					#				  "sphingoid base - CH2O - ADDCUT"])

					FRAGMENTS.append([MW("C11H19NO9") + PROTON, 7, "N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9") - H2O + PROTON, 5, "deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9") - 2 * H2O + PROTON, 10, "deoxy-N-Ac-Neuraminic Acid - water "])

					FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 6, "deoxy-N-Ac-glucosamine"])
					FRAGMENTS.append([MW("C8H15NO6") - H2O + ADDUCT[adduct], 8, "deoxy-N-Ac-glucosamine + adduct"])
					FRAGMENTS.append([MW("C8H15NO6") - 2*H2O + ADDUCT[adduct], 17, "deoxy-N-Ac-glucosamine - water + adduct"])

					FRAGMENTS.append([MW("C6H12O6") + MW("C11H19NO9") - 2 * H2O + PROTON, 35, "deoxygalactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3 * H2O + PROTON, 3,
									  "deoxy-N-Ac-galcatosamine + deoxy-galactose + deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C6H12O6") + MW("C8H15NO6") + MW("C11H19NO9") - 3 * H2O + ADDUCT[adduct], 11,
									  "deoxy-N-Ac-galcatosamine + deoxy-galactose + deoxy-N-Ac-Neuraminic Acid + ADDUCT"])

					FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + PROTON, 2, "deoxy-N-Ac-glucosamine + deoxygalactose"])
					FRAGMENTS.append([MW("C8H15NO6") + MW("C6H12O6") - 2 * H2O + ADDUCT[adduct], 105, "deoxy-N-Ac-glucosamine + deoxygalactose + ADDUCT"])



			
			elif adduct in NEG_ADDUCTS:  
				if adduct == "[M+Cl]-" or adduct == "[M+FA]-":
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON, 1000, "NL adduct"] )
				if adduct == "[M-H]-":
					FRAGMENTS.append([PREC, 15, "precursor"])
					#FRAGMENTS.append( [PREC - H2O, 5, "NL water"] )
					FRAGMENTS.append( [PREC - MW("CO2"), 7, "NL CO2"] )
					FRAGMENTS.append( [MW("C11H19NO9")-PROTON, 3, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 1000, "NL  deoxy-N-Ac-Neuraminic Acid"] )
					#FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("CO2")+H2O, 4, "NL deoxy-N-Ac-Neuraminic Acid + CO2"] )
					FRAGMENTS.append( [PREC-2*MW("C11H19NO9")+2*H2O, 3, "NL 2x deoxy-N-Ac-Neuraminic Acid"] )

					#FRAGMENTS.append( [PREC-MW("C11H19NO9"), 5, "NL N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON, 259, "N-acyl-ethanolamine - H2O - H"])
					#FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON, 2, "Fatty Amide - NH3 - H"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H10O5") + 3 * H2O, 10, "Ceramide"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 2 * H2O, 2, "Glucosyl ceramide - water"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 3* H2O, 9, "Glucosyl ceramide"])
					FRAGMENTS.append([PREC - 2 * MW("C11H19NO9") - MW("C8H15NO6") -  MW("C6H10O5") + 3 * H2O, 41, "Lactosyl ceramide"])

				if adduct == "[M-2H]2-":
					SINGLE = MASS - PROTON
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append([PREC - H2O / 2, 13, "NL water 2-"])
					FRAGMENTS.append([PREC - MW("CO2") / 2, 55, "NL CO2 2-"])
					FRAGMENTS.append([MW("C11H19NO9") - PROTON, 2, "N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9") - H2O - PROTON, 343, "deoy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([SINGLE - MW("C11H19NO9") + H2O, 133, "NL  deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") + 2 * H2O, 6, "NL 2x deoxy-N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([SINGLE - MW("C11H19NO9"), 4, "NL N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([NL(self.chains[1]) - H2O + MW("C2H7NO") - H2O - PROTON, 9, "N-acyl-ethanolamine - H2O - H"])

					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 3 * MW("C6H10O5") + 3 * H2O, 13, "Ceramide"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 2 * H2O, 3, "Glucosyl ceramide - water"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - 2 * MW("C6H10O5") + 3 * H2O, 7, "Glucosyl ceramide"])
					FRAGMENTS.append([SINGLE - 2 * MW("C11H19NO9") - MW("C8H15NO6") - MW("C6H10O5") + 3 * H2O, 24, "Lactosyl ceramide"])

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




x  = AcGD1a("AcGD1a",[[18,1,1], [18,0,0]],  adduct="[M+H]+")
print x.printNist()



x  = AcGD1a("AcGD1a",[[18,1,1], [18,0,0]],  adduct="[M+2H]2+")
print x.printNist()


x  = AcGD1a("AcGD1a",[[18,1,1], [18,0,0]],  adduct="[M+Na]+")
print x.printNist()
