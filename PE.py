from lipidRules import *

class PE(GPL):




		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
		
		
				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 1, "precursor"])
				
					FRAGMENTS.append( [PREC-MW("C2H8NO4P"), 1000, "NL ethanolamine-P"] )

# check assymetry between sn1 and sn2 in reference spectra

					i = 0
					chain = str(i+1)
					FRAGMENTS.append( [PREC - NL(self.chains[i]), 1, "NL sn" + chain ])
					FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 40, "sn" + chain + " acylium ion - water"])
					FRAGMENTS.append( [NL(self.chains[i])+PROTON, 13, "sn" + chain + " acylium ion"])
					#below is probably not a real peak, but added for DIMS ease
					FRAGMENTS.append([PREC - MW("C2H5N") - NL(self.chains[i]), 5, "NL sn" + chain + " + ethanolamine"])

					i = 1
					chain = str(i+1)
					FRAGMENTS.append( [PREC - NL(self.chains[i]), 4, "NL sn" + chain ])
					FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 26, "sn" + chain + " acylium ion - water"])
					FRAGMENTS.append( [PREC-MW("C2H5N")-NL(self.chains[i]), 2, "NL sn" + chain + " + ethanolamine"] )
					# below is probably not a real peak, but added for DIMS ease
					FRAGMENTS.append([NL(self.chains[i]) + PROTON, 5, "sn" + chain + " acylium ion"])

				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 374, "precursor"])
					FRAGMENTS.append( [MW("PO4H3")+ADDUCT[adduct], 900, "phosphoric acid + adduct"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct], 1000, "C2H8NO4P + adduct"] )
					FRAGMENTS.append( [PREC-MW("C2H5N"), 225, "NL ethanolamine"] )
					FRAGMENTS.append( [PREC-MW("C2H8NO4P"), 645, "NL ethanolamine-P"] )
					FRAGMENTS.append( [PREC-MW("C2H8NO4P")+H2O, 27, "NL deoxyethanolamine-P"] )
					FRAGMENTS.append( [PREC-MW("C2H8NO4P")-ADDUCT[adduct]+PROTON, 204, "NL ethanolamine-P + adduct"] )
	  
					i = 0
					chain = str(i+1)
					FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O-H2O, 2, "sn" + chain + " acylium ion - water"])
					FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 20, "sn" + chain + " acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]), 7, "NL sn" + chain ])
					FRAGMENTS.append( [PREC-MW("C2H5N")-NL(self.chains[i]), 36, "NL sn" + chain + " + ethanolamine"] ) 

					i =1
					chain = str(i+1)
					FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O-H2O, 9, "sn" + chain + " acylium ion - water"])
					FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 22, "sn" + chain + " acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]), 2, "NL sn" + chain ])
					FRAGMENTS.append( [PREC-MW("C2H5N")-NL(self.chains[i]), 23, "NL sn" + chain + " + ethanolamine"] ) 


		
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 111, "precursor"])
				FRAGMENTS.append( [MW("C3H9O6P")-H2O-PROTON, 20,  "glycerol-P - water"] )
				FRAGMENTS.append( [MW("H3PO4")-PROTON,  5,  "phosphate"] ) 
				FRAGMENTS.append( [MW("HPO3")-PROTON,  49,  "phosphite"] ) 
				FRAGMENTS.append( [MW("C5H12NO5P")-PROTON, 27, "ethanaolamine-P + deoxyglycerol"] )

				FRAGMENTS.append( [MW("C2H8NO4P")-H2O-PROTON,  9,  "deoxyethanolamine-P"] )
				FRAGMENTS.append( [MW("C2H8NO4P")-PROTON,  43,  "ethanolamine-P"] )

				i = 0
				chain = str(i+1)
				FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 8, "NL sn" + chain + "ketene"])
				FRAGMENTS.append( [PREC - NL(self.chains[i]), 5, "NL sn" + chain])
				if self.chains[i][1] > 3:
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 969, "sn" + chain])
					FRAGMENTS.append( [NL(self.chains[i])-MW("CO2")-PROTON, 93, "sn" + chain + "decarboxylation"])
				else :
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn" + chain])

				i = 1
				chain = str(i+1)
				FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 45, "NL sn" + chain + "ketene"])
				FRAGMENTS.append( [PREC - NL(self.chains[i]), 10, "NL sn" + chain])
				if self.chains[i][1] > 3:
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 969, "sn" + chain])
					FRAGMENTS.append( [NL(self.chains[i])-MW("CO2")-PROTON, 93, "sn" + chain + "decarboxylation"])
				else: 
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn" + chain])



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




x  = PE("PE",[[16,0,0], [18,1,0]],  adduct="[M-H]-")
print x.printNist()

