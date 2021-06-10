from lipidRules import *

class APCS(singleAcyl):
#based upon Sidu et al JLR 2019

		neg_adduct_set = ["[M-H]-"]
		pos_adduct_set = ["[M+H]+"]

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
	


				if adduct == "[M+H]+":  
					FRAGMENTS.append([PREC, 10, "precursor"])
					FRAGMENTS.append([PREC - MW("C5H14NO4P"), 5, "NL Phosphocholine"])
					FRAGMENTS.append( [MW("C5H13NO")+PROTON-H2O, 80, "Dehydrocholine"] )
					FRAGMENTS.append( [MW("C2H5O4P")+PROTON, 40, "C2H6O4P"] )
					FRAGMENTS.append( [MW("C5H14NO4P")+PROTON, 1000, "Phosphocholine"] )

					FRAGMENTS.append( [PREC - NL(self.chains[0])+H2O-MW("CO2H2"), 10, "NL Acyl + Formic Acid"] )
					FRAGMENTS.append( [PREC - NL(self.chains[0])+H2O, 6, "NL Acyl"] )

				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 5, "precursor"])
					FRAGMENTS.append([PREC - MW("C5H14NO4P"), 280, "NL Phosphocholine"])
					FRAGMENTS.append( [MW("C5H13NO")+PROTON-H2O, 160, "Dehydrocholine"] )
					FRAGMENTS.append( [MW("C5H14NO4P")+ADDUCT[adduct], 1000, "Phosphocholine + Addcut"] )
					FRAGMENTS.append( [MW("C2H5O4P")+ADDUCT[adduct], 700, "C2H6O4P + Adduct"] )
					FRAGMENTS.append([PREC - MW("C3H9N"), 75, "NL trimethylamine"])



			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 80, "precursor"])
				FRAGMENTS.append([PREC-MW("CO2")-MW("C5H14NO4P"), 10, "NL CO2 + phosphocholine"])
				FRAGMENTS.append([PREC-MW("C5H14NO4P"), 10, "NL phosphocholine"])
				FRAGMENTS.append([MW("HPO3") - PROTON, 200, "Phosphite"])
				FRAGMENTS.append([MW("H3PO4") - PROTON, 15, "Phosphate"])
				FRAGMENTS.append([MW("C5H15NO4P") - PROTON - MW("CH3"), 220, "PC - CH3"])
				FRAGMENTS.append([NL(self.chains[0])-H2O - PROTON, 25, "ketene"])
				FRAGMENTS.append([MW("C2H5O4P") - PROTON, 1000, "C2H6O4P"])

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




x  = APCS("APCS",[[16,0,0]],  adduct="[M+Na]+")
print x.printNist()

