from lipidRules import *

class DG(GPL):

		pos_adduct_set = ["[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+Li]+"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 1000, "precursor"])


				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 1000, "precursor"])

					for i in range(0,len(self.chains)):
						chain = str(i+1)			
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O+MW("C3H6O2"), 14, "sn"+chain+" + deoxyglycerol"] )
 						FRAGMENTS.append( [NL(self.chains[i])+ADDUCT[adduct]-H2O+MW("C3H6O2"), 5, "sn"+chain+" + deoxyglycerol + adduct"] )


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




x  = DG("DG",[[16,0,0], [20,1,0]],  adduct="[M+NH4]+")
print x.printNist()
