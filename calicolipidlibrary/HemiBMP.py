from lipidRules import *

class HemiBMP(LysoCardioLipin):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]


			
			if adduct in NEG_ADDUCTS:
				FRAGMENTS.append([PREC, 447, "precursor"])

				FRAGMENTS.append( [MW("C3H7O5P")-PROTON, 41, "C3H6O5P-"] )



				for i in range(0,len(self.chains)):
					chain = str(i + 1);		
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 333, "sn"+chain+" RCO2-"])
					FRAGMENTS.append( [PREC-NL(self.chains[i]), 3, "NL (sn"+chain+" + water)"])
					FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O, 4, "NL sn"+chain])

					FRAGMENTS.append( [NL(self.chains[i])+MW("C3H7O5P")-PROTON, 1, "sn"+chain+" LPA"])
					FRAGMENTS.append( [NL(self.chains[i])+MW("C3H7O5P")-H2O-PROTON, 8, "sn"+chain+" LPA - water"])


								

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





x  = HemiBMP("HemiBMP",[[18,0,0], [22,6,0], [2,0,0]],  adduct="[M-H]-")
print x.printNist()

