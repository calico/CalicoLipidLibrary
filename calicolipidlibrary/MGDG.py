from lipidRules import *

class MGDG(GPL):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+Na]+" or adduct == "[M+K]+":
					FRAGMENTS.append([PREC, 1000, "precursor"])
	
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H10O5")+PROTON, 1000, "NL hexose + adduct"] )
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H10O5")-H2O+PROTON, 1000, "NL hexose + adduct + water"] )
					FRAGMENTS.append( [MW("C9H16O6")+ADDUCT[adduct], 100, "hexosylglycerol - 2x water"] )

				for i in range(0,len(self.chains)):
						chain = str(i+1)
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O+MW("C3H6O2"), 1000, "sn"+chain+" + glycerol - water"] )
						FRAGMENTS.append( [PREC - NL(self.chains[i]), 1000, "NL sn"+chain] )
				
				if adduct == "[M+NH4]+":  
					FRAGMENTS.append([PREC, 1000, "precursor"])

					FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H10O5")+PROTON, 1000, "NL hexose + adduct"] )
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H10O5")-H2O+PROTON, 1000, "NL hexose + water + adduct"] )

					for i in range(0,len(self.chains)):
						chain = str(i+1)
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O+MW("C3H6O2"), 1000, "sn"+chain+" + glycerol- water"] )
 
			if adduct in NEG_ADDUCTS:
				FRAGMENTS.append([PREC, 1000, "precursor"])
				for i in range(0,len(self.chains)):
					chain = str(i+1)
					FRAGMENTS.append([NL(self.chains[i])-PROTON, 1000, "sn" + chain])	
				


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




x  = MGDG("MGDG",[[18,3,0], [16,3,0]],  adduct="[M+NH4]+")
print x.printNist()


