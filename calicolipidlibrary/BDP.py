from lipidRules import *

class BDP(CardioLipin):

		neg_adduct_set = ["[M-H]-"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]


			if adduct in POS_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])


			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 458, "precursor"])

				FRAGMENTS.append( [MW("C3H7O5P")-PROTON, 4, "C3H6O5P-"] )

				for i in range(0,len(self.chains)):
					chain = str(i + 1);		
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 250, "sn"+chain+" RCO2-"])
					FRAGMENTS.append( [PREC-NL(self.chains[i]), 2, "NL (sn"+chain+")"])
					FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O, 3, "NL sn"+chain+" ketene"])

					FRAGMENTS.append( [NL(self.chains[i])+MW("C3H7O5P")-PROTON, 1, "sn"+chain+" LysoPA"])
					FRAGMENTS.append( [NL(self.chains[i])+MW("C3H7O5P")-H2O-PROTON, 2, "sn"+chain+" LysoPA - water"])


								

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





# x  = BDP("BDP",[[18,1,0], [18,1,0], [18,1,0],[18,1,0]],  adduct="[M-H]-")
# print x.printNist()

# calicolipidlibrary.print_spectrum("BDP", [[18,1,0], [18,0,0], [16,1,0],[16,0,0]], "[M-H]-")