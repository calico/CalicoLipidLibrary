from lipidRules import *

class MMPE(GPL):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
	

		
				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 2, "precursor"])
					FRAGMENTS.append([MW("C3H10NO4P") + ADDUCT[adduct], 35, "MMPE phosphate"])
					FRAGMENTS.append([MW("C3H7N") + PROTON, 350, "MMPE - water"])
					FRAGMENTS.append( [PREC-MW("C3H10NO4P"), 1000, "NL MMPE + phosphate"] )
	
				
					i=0
					chain = str(i+1)
					#FRAGMENTS.append( [PREC - NL(self.chains[i])-H2O, 1000, "NL sn" + chain + " ketene"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]), 1, "NL sn" + chain ])
					FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 40, "sn" + chain + " acylium ion "])
					FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O-H2O, 15, "sn" + chain + " acylium ion - water"])
					FRAGMENTS.append( [NL(self.chains[i])+PROTON+MW("C3H4O"), 10, "sn" + chain + " acylium ion + deoxyglycerol"])

					i=1
					chain = str(i + 1)
					#FRAGMENTS.append([PREC - NL(self.chains[i]) - H2O, 1000, "NL sn" + chain + " ketene"])
					FRAGMENTS.append([PREC - NL(self.chains[i]), 1, "NL sn" + chain])
					FRAGMENTS.append([NL(self.chains[i]) + PROTON - H2O, 85, "sn" + chain + " acylium ion "])
					FRAGMENTS.append([NL(self.chains[i]) + PROTON - H2O - H2O, 32, "sn" + chain + " acylium ion - water"])
					FRAGMENTS.append([NL(self.chains[i]) + PROTON + MW("C3H4O"), 10, "sn" + chain + " acylium ion + deoxyglycerol"])


				if adduct == "[M+Na]+" or adduct == "[M+NH4]+" or adduct == "[M+K]+":
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append([MW("C3H10NO4P") + ADDUCT[adduct], 1000, "MMPE phosphate"])
					FRAGMENTS.append([MW("C3H7N") + PROTON, 1000, "MMPE - water"])
					FRAGMENTS.append( [MW("PO4H3")+ADDUCT[adduct], 1000, "phosphoric acid + adduct"] )
					FRAGMENTS.append( [PREC-MW("C3H7N"), 1000, "NL MMPE"] )
					FRAGMENTS.append( [PREC-(MW("C3H10NO4P")), 1000, "NL MMPE + phosphate"] )
					FRAGMENTS.append( [PREC-MW("C3H10NO4P")-ADDUCT[self.adduct]+PROTON, 1000, "NL MMPE + phosphate + ADDUCT"] )
				
	  
					for i in range(0,len(self.chains)):
						chain = str(i+1)
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O-H2O, 1000, "sn" + chain + " acylium ion - water"])
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 1000, "sn" + chain + " acylium ion"])
						FRAGMENTS.append( [PREC - NL(self.chains[i]), 1000, "NL sn" + chain + "" ])
						FRAGMENTS.append( [PREC-MW("C3H7N")-NL(self.chains[i]), 1000, "NL sn" + chain + " + MMPE"] ) 

		
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])



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




# x  = MMPE("MMPE",[[16,1,0], [16,1,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("MMPE",[[16,1,0], [16,1,0]], "[M+H]+")