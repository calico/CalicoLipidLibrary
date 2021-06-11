from lipidRules import *

class PG(GPL):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]
			FRAGMENTS.append([PREC, 1000, "precursor"])


			if adduct in POS_ADDUCTS:  


				if adduct == "[M+H]+":
				
					FRAGMENTS.append( [PREC - MW("C3H9O6P"), 1000, "NL glycerol-P"])

					for i in range(0,len(self.chains)):
						chain = str(i + 1);
				
						FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 1000, "RCO+ sn"+chain])
						FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O+H2O, 1000, "NL sn"+chain+"-H2O"])


				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([MW("C3H9O6P")+ADDUCT[adduct], 1000, "glycerol-P + adduct"])
					FRAGMENTS.append([MW("C3H9O6P")-H2O+ADDUCT[adduct], 43, "deoxyglycerol-P + adduct"])

					
					FRAGMENTS.append( [PREC - MW("C3H9O6P"), 43, "NL glycerol-P"])
					FRAGMENTS.append( [PREC - MW("C3H9O6P")-ADDUCT[adduct]+PROTON, 79, "NL glycerol-P + Adduct"])

					i = 0
					chain = str(i + 1)
					FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 5, "sn"+chain+" acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 1, "NL sn"+chain+" ketene"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]), 6, "NL sn"+chain+""])


					i = 1
					chain = str(i + 1)
					FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 7, "sn"+chain+" acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 3, "NL sn"+chain+" ketene"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]), 3, "NL sn"+chain+""])


		
			elif adduct in NEG_ADDUCTS:  

				FRAGMENTS.append( [MW("H3PO4")-PROTON,  5,  "phosphate"] ) 
				FRAGMENTS.append( [MW("HPO3")-PROTON,  40,  "phosphite"] ) 
				FRAGMENTS.append( [MW("C3H7O5P")-PROTON, 116, "deoxyglycerol-P-"] )
				FRAGMENTS.append( [MW("C3H9O6P")-PROTON, 34, "glycerol-P"] )
				FRAGMENTS.append( [MW("C6H13O7P")-PROTON, 13, "NL sn1 + sn2"] )
				FRAGMENTS.append( [MW("C6H11O6P")-PROTON, 6, "NL sn1 + sn2 + water"] )

				i = 0
				chain = str(i + 1)
				FRAGMENTS.append( [NL(self.chains[i])-PROTON, 530, "sn1"])
				FRAGMENTS.append( [PREC-NL(self.chains[i]), 7, "nL (sn1 + water)"])
				FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O, 5, "nL sn1 "])
				FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O-MW("C3H8O3"), 8, "nL (sn1 + glycerol)"])

				i = 1
				chain = str(i + 1)	
				FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn2"])
				FRAGMENTS.append( [PREC-NL(self.chains[i]), 13, "nL (sn2 + water)"])
				FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O, 18, "nL sn2 "])
				FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O-MW("C3H8O3"), 22, "nL (sn2 + glycerol)"])

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




x  = PG("PG",[[16,0,0], [18,1,0]],  adduct="[M+Na]+")
print x.printNist()

