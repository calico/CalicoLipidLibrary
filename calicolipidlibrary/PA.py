from lipidRules import *

class PA(GPL):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

		
		
			if adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 291, "precursor"])

				FRAGMENTS.append( [MW("H3PO4")-PROTON,  23,  "phosphate"] ) 
				FRAGMENTS.append( [MW("C3H9O6P")-PROTON, 6, "glycerol phosphate"] )
				FRAGMENTS.append( [MW("C3H7O5P")-PROTON, 302, "glycerol phosphate - water"] )
				FRAGMENTS.append( [MW("C3H5O4P")-PROTON, 6, "glycerol phosphate - 2x water"] )
	
				for i in range(0,len(self.chains)):
					chain = str(i + 1);
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 500, "sn"+chain])
					FRAGMENTS.append( [PREC-NL(self.chains[i]), 138, "nL sn"+chain])
					FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O, 52, "nL sn"+chain+" ketene"])

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



# x  = PA("PA",[[16,0,0], [18,1,0]],  adduct="[M-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("PA",[[16,0,0], [18,1,0]], "[M-H]-")