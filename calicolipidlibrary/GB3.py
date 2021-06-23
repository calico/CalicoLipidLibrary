from lipidRules import *

class GB3(SphingoLipid):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":    
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append( [PREC-H2O, 1000, "NL H2O"] )
					FRAGMENTS.append( [PREC-MW("C6H12O6"), 1000, "NL hexose"])

					FRAGMENTS.append( [PREC-MW("C12H22O11"), 1000, "NL 2x galactose"])
					FRAGMENTS.append( [PREC-MW("C12H22O11")-H2O, 1000, "NL 2x galactose + water"])
					FRAGMENTS.append( [PREC-MW("C18H32O16"), 1000, "Ceramide"])
					FRAGMENTS.append( [PREC-MW("C18H32O16")-H2O, 1000, "Ceramide - water"])
					FRAGMENTS.append( [PREC - MW("C18H32O16")-NL(self.chains[1])+H2O, 1000, "sphingobase"])
					FRAGMENTS.append( [PREC - MW("C18H32O16")-NL(self.chains[1]), 1000, "sphingobase - water"])
					FRAGMENTS.append( [PREC - MW("C18H32O16")-NL(self.chains[1])+H2O-MW("CH2O"), 1000, "sphingobase - CH2O"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 1000, "fatty amide"])


				if adduct == "[M+Na]+" or adduct == "[M+K]+" or adduct == ["M+NH4]+"]:
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append( [PREC-MW("C6H12O6")+H2O, 1000, "NL deoxygalactose"])
					FRAGMENTS.append( [PREC - NL(self.chains[1]) - MW("C18H32O16")  - ADDUCT[adduct] + PROTON, 10, "N'': NL 3x hexose + acyl ketene + water + ADDUCT"] )
					FRAGMENTS.append( [MW("C18H32O16")+ADDUCT[adduct], 10, "2x galactose + glucose + Adduct"])
					FRAGMENTS.append( [MW("C18H32O16")-H2O+ADDUCT[adduct], 10, "2x galactose + deoxyglucose + Adduct"])


			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])

				if adduct == "[M+Cl]-" or adduct == "[M+FA]-":
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON, 1000, "NL adduct"] )
				if adduct == "[M-H]-":
					FRAGMENTS.append( [PREC-MW("C6H10O5"), 1000, " NL deoxyhexose"])
				
					FRAGMENTS.append( [PREC-MW("C6H12O6"), 1000, "NL hexose"])

					FRAGMENTS.append( [PREC-MW("C12H20O10"), 1000, "NL deoxylactose"])
					FRAGMENTS.append( [PREC-MW("C12H22O11"),1000 , "NL lactose"])

					FRAGMENTS.append([PREC - MW("C18H32O16"), 1000, "NL 2x galactose + deoxyglucose"])
					FRAGMENTS.append([PREC - MW("C18H32O16"), 1000, "NL 2x galactose + deoxyglucose"])

					FRAGMENTS.append( [MW("C12H22O11")-H2O-PROTON, 1000, "deoxylactose"])
					FRAGMENTS.append( [MW("C18H32O16")-H2O-PROTON, 1000, "2x galactose + deoxyglucose"])
					FRAGMENTS.append( [MW("C18H32O16")-PROTON, 1000, "2x galactose + glucose"])


					FRAGMENTS.append( [MW("C6H12O6")-PROTON, 1000, "hexose"])

					#FRAGMENTS.append( [MW("C5H6O3")-PROTON, 167, "hexose fragment at MZ 113"])
					#FRAGMENTS.append( [MW("C4H8O4")-PROTON, 61, "hexose fragment at MZ 119"])
					#FRAGMENTS.append( [MW("C5H8O4")-PROTON, 38, "hexose fragment at MZ 131"])
					#FRAGMENTS.append( [MW("C6H8O4")-PROTON, 79, "hexose fragment at MZ 143"])
					#FRAGMENTS.append( [MW("C5H10O5")-PROTON, 13, "hexose fragment at MZ 149"])
					#FRAGMENTS.append( [MW("C6H10O5")-PROTON, 83, "hexose fragment at MZ 161"])
					
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")-PROTON, 1000, "sn2 fatty amide"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+MW("C2H2")-PROTON, 1000, "sn2 fatty amide + C2H2"])
					FRAGMENTS.append( [NL(self.chains[1])-PROTON, 1000, "sn2 fatty acid"])



				if (self.chains[1][2] == 1): #if 2-hydroxy FA
					FRAGMENTS.append( [ NL(self.chains[1]) - H2O -  MW("CO") - PROTON, 1000, "Fatty Acyl - CH2O2"])  #fragment intensity not from standard
					FRAGMENTS.append( [ NL(self.chains[1]) - H2O + MW("H2")- PROTON, 1000, "Fatty Acyl +H2"])	#fragment intensity not from standard

				if (self.chains[0][2] > 0): #if not a m LCB
					FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("C18H32O16") - MW("C2H5N") - PROTON, 1000, "LCB-CO2H4-CNH3"])


				if (self.chains[0][2] == 2): #it a t LCB
					FRAGMENTS.append( [PREC - ADDUCT[adduct] - NL(self.chains[1]) +H2O - MW("C18H32O16") - MW("C3H7NO") - PROTON, 1000, "LCB-CO2H4-CNH3 - (CH2O if t) "])


				

			
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



# x  = GB3("GB3",[[18,1,1], [16,0,0]],  adduct="[M+H]+")
# print x.printNist()


