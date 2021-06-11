from lipidRules import *

class HexCer(SphingoLipid):

		pos_adduct_set = ["[M+H]+", "[M+Na]+", "[M+K]+"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 4, "precursor"])

					if self.chains[0][1] != 0:
						FRAGMENTS.append( [PREC-H2O, 73, "NL H2O"] )
						FRAGMENTS.append( [PREC-MW("C6H12O6"), 117, "NL hexose"])
						FRAGMENTS.append( [PREC-MW("C6H12O6")-H2O, 65, "NL hexose  + water"])
						FRAGMENTS.append( [PREC - MW("C6H12O6")-NL(self.chains[1])+H2O, 95, "NL hexose + acyl ketene"])
						FRAGMENTS.append( [PREC - MW("C6H12O6")-NL(self.chains[1]), 1000, "NL hexose + acyl ketene + water"])
						FRAGMENTS.append( [PREC - MW("C6H12O6")-NL(self.chains[1])+H2O-MW("CH2O"), 93, "NL hexose + acyl ketene + water + CH2O)"])
						FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 21, "fatty amide"])
						if (self.chains[1][2] > 0):
							 FRAGMENTS.append( [NL(self.chains[1])-H2O-H2O+MW("NH3")+PROTON, 20, "fatty amide - water"]) #guess, not from spectra
					else: #for dihydrosphingosine backbones
						FRAGMENTS.append([PREC - H2O, 2, "NL H2O"])
						FRAGMENTS.append([PREC - MW("C6H12O6") + H2O, 337, "NL deoxyhexose"])
						FRAGMENTS.append([PREC - MW("C6H12O6"), 1000, "NL hexose"])
						FRAGMENTS.append([PREC - MW("C6H12O6") - H2O, 59, "NL hexose  + water"])
						FRAGMENTS.append([PREC - MW("C6H12O6") - NL(self.chains[1]) + H2O, 391, "NL deoxyhexose + acyl ketene"])
						FRAGMENTS.append([PREC - MW("C6H12O6") - NL(self.chains[1]) + 2*H2O, 81, "NL double-deoxyhexose + acyl ketene"])

						FRAGMENTS.append([PREC - MW("C6H12O6") - NL(self.chains[1]), 232, "NL hexose + acyl ketene"])
						FRAGMENTS.append([PREC - MW("C6H12O6") - NL(self.chains[1]) + H2O - MW("CH2O"), 27, "NL hexose + acyl ketene + water + CH2O)"])
						FRAGMENTS.append([NL(self.chains[1]) - H2O + MW("NH3") + PROTON, 259, "fatty amide"])
						if (self.chains[1][2] > 0):
							FRAGMENTS.append([NL(self.chains[1]) - H2O - H2O + MW("NH3") + PROTON, 50, "fatty amide - water"]) # guess, not from spectra

				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O+PROTON, 1000, "sn2 acylium ion"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O-H2O+PROTON, 1000, "sn2 acylium ion - H2O"])
					FRAGMENTS.append([PREC - MW("C6H12O6"), 1000, "NL hexose"])
					FRAGMENTS.append([PREC - MW("C6H12O6")+H2O, 1000, "NL deoxyhexose"])
					FRAGMENTS.append([PREC - MW("C6H12O6")-MW("COH2"), 1000, "NL hexose + H2CO"])
					FRAGMENTS.append([PREC - MW("C6H12O6") - NL(self.chains[1]), 1000, "NL hexose + acyl ketene + water"])
					FRAGMENTS.append([MW("C6H12O6") + ADDUCT[adduct], 1000, "hexose + adduct"])


			elif adduct in NEG_ADDUCTS:

				if adduct == "[M-H]-":
					FRAGMENTS.append([PREC, 452, "precursor"])

				elif adduct == "[M+Cl]-" or adduct == "[M+FA]-" or adduct == "[M+37Cl]-":
					FRAGMENTS.append([PREC, 1, "precursor"])
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON, 452, "pre - adduct"] )

				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H10O5")-PROTON, 766, "NL deoxyhexose + adduct"])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H12O6")-PROTON, 21, "NL hexose + adduct"])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C7H12O6")-PROTON, 41, "NL C7H12O6 + adduct"])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C7H14O6")-PROTON, 22, "NL C7H14O6 + adduct"])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C7H12O6")-H2O-PROTON, 45, "NL C7H12O6 + H2O + adduct"])

				FRAGMENTS.append( [MW("C6H12O6")-PROTON, 137, "hexose"])
				FRAGMENTS.append( [MW("C3H6O3")-PROTON, 911, "hexose fragment at MZ 89"])
				FRAGMENTS.append( [MW("C4H6O3")-PROTON, 769, "hexose fragment at MZ 101"])
				FRAGMENTS.append( [MW("C5H6O3")-PROTON, 379, "hexose fragment at MZ 113"])
				FRAGMENTS.append( [MW("C4H8O4")-PROTON, 424, "hexose fragment at MZ 119"])
				FRAGMENTS.append( [MW("C5H8O4")-PROTON, 55, "hexose fragment at MZ 131"])
				FRAGMENTS.append( [MW("C6H8O4")-PROTON, 77, "hexose fragment at MZ 143"])
				FRAGMENTS.append( [MW("C5H10O5")-PROTON, 33, "hexose fragment at MZ 149"])
				FRAGMENTS.append( [MW("C6H10O5")-PROTON, 71, "hexose fragment at MZ 161"])

				FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")-PROTON, 84, "sn2 fatty amide"])
				FRAGMENTS.append( [NL(self.chains[1])-PROTON, 128, "sn2 fatty acid"])
				FRAGMENTS.append([NL(self.chains[1]) - H2O - PROTON, 71, "acyl ketene"])

				FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+MW("C2H2")-PROTON, 309, "sn2 fatty amide + C2H2"])
				FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+MW("C2H2O")-PROTON, 65, "sn2 fatty amide + C2H2O"])

				FRAGMENTS.append([PREC - ADDUCT[adduct] - NL(self.chains[1]) + H2O - MW("C6H12O6 ") - MW("NH3") - PROTON, 62, "LCB - NH3 "])

				if (self.chains[1][2] == 1):
					FRAGMENTS.append( [ NL(self.chains[1]) - H2O -  MW("CO") - PROTON, 10, "acyl ketene - CH2O2"]) #intensity a guess, not from spectra
					FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("H2")- PROTON, 10, "acyl ketene + H2"]) #intensity a guess, not from spectra

#					FRAGMENTS.append( [ NL(self.chains[1]) - H2O - PROTON + MW("H2") - MW("O"), 1000, "Fatty Amide - O + 2H - NH3 - H"])


				if (self.chains[0][2] > 0): 
					#FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("C6H12O6 ")  -MW("NC2H3") - PROTON, 1000, "LCB - NC2H3 "])
					FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("C6H12O6 ") - PROTON, 3, "LCB  "])
					FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("C6H12O6") - MW("C2H5N") - PROTON, 95, "LCB-CO2H4-CNH3"])


					
				if (self.chains[0][2] == 2): 
					FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("C6H12O6") - MW("C3H7NO") - PROTON, 10, "LCB-CO2H4-CNH3 - (CH2O if t) "]) #intensity a guess, not grom spectra


				

			
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




x  = HexCer("HexCer",[[18,1,2], [18,2,0]],  adduct="[M+H]+")
print x.printNist()

