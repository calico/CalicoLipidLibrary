from lipidRules import *

class LysoSM(LysoSphingoLipid):

		neg_adduct_set = GPL.neg_adduct_set + ["[M-CH3]-"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
	


				if adduct == "[M+H]+":  
					FRAGMENTS.append([PREC, 92, "precursor"])
			
					FRAGMENTS.append( [MW("C5H13NO")+PROTON, 18, "Choline"] )
					FRAGMENTS.append( [MW("C5H13NO")+PROTON-H2O, 35, "Dehydrocholine"] )
					FRAGMENTS.append( [MW("C5H14NO4P")+PROTON, 1000, "Phosphocholine"] )
			
					FRAGMENTS.append( [PREC-H2O, 13, "NL water"] )
					FRAGMENTS.append( [MW("C2H5O4P")+PROTON, 15, "C2H6O4P"] )
					FRAGMENTS.append( [MW("C3H9N")+ADDUCT[adduct], 5, "Trimethylamine"] )

					if (self.chains[0][1] == 0):
						FRAGMENTS.append( [PREC-MW("C5H14NO4P")-H2O, 11, "NL Phosphocholine + water"] )   	
						FRAGMENTS.append( [PREC-MW("C5H14NO4P"), 6, "NL Phosphocholine"] )   	
					else:
						FRAGMENTS.append( [PREC-MW("C5H14NO4P")-H2O, 20, "NL Phosphocholine + water"] )   	
				
			

				if adduct == "[M+Na]+" or adduct == "[M+K]+" or adduct == "[M+NH4]+":
					FRAGMENTS.append([PREC, 1000, "precursor"])

					FRAGMENTS.append( [MW("C2H5O4P")+ADDUCT[adduct],920, "C2H5O4P + Adduct"] )
					FRAGMENTS.append( [PREC - MW("C3H9N"), 375, "NL Trimethylamine"] )
					FRAGMENTS.append( [PREC - MW("C5H14NO4P"), 639, "NL PhosphoCholine"] )

					FRAGMENTS.append( [MW("C5H14NO"), 108, "Choline"] )
					FRAGMENTS.append( [MW("C5H14NO")-H2O, 201, "Dehydrocholine"] )  



		
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 92, "precursor"])
		
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON-MW("CH3"), 1000, "pre-adduct-CH3"] )
		
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




x  = LysoSM("LysoSM",[[16,0,0]],  adduct="[M+H]+")
print x.printNist()



x  = LysoSM("LysoSM",[[16,0,0]],  adduct="[M+Na]+")
print x.printNist()

