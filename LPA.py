from lipidRules import *

class LPA(LysoGPL):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			
			
			if adduct in NEG_ADDUCTS:

				if self.chains[0][0] != 0:
					i = 0
				else:
					i = 1
				chain = str(i + 1)

				FRAGMENTS.append([PREC, 54, "precursor"])
				
				FRAGMENTS.append( [MW("H3PO4")-PROTON,  20,  "phosphate"] ) 
				FRAGMENTS.append( [MW("HPO3")-PROTON,  223,  "phosphite"] ) 
				FRAGMENTS.append( [MW("C3H5O5P")-PROTON, 6, "glycerol phosphate - water - H2"] )
				FRAGMENTS.append( [MW("C3H7O5P")-PROTON, 1000, "glycerol phosphate - water"] )
				FRAGMENTS.append( [MW("C3H9O6P")-PROTON, 32, "glycerol phosphate"] )
		
				FRAGMENTS.append( [NL(self.chains[i])-PROTON, 17, "sn"+ chain])

								


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




x  = LPA("LPA",[[16,0,0],[0,0,0]],  adduct="[M-H]-")
print x.printNist()
