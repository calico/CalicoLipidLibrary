from lipidRules import *

class LCB_P(LysoSphingoLipid):
		neg_adduct_set = ["[M-H]-"]
		pos_adduct_set = ["[M+H]+"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":  
					FRAGMENTS.append([PREC, 3, "precursor"])
					FRAGMENTS.append( [PREC-H2O, 13, "NL water"] )
					FRAGMENTS.append( [PREC-MW("PO4H3"), 6, "NL phosphoric acid"] )
					
					FRAGMENTS.append( [PREC-MW("PO4H3")-H2O, 1000, "NL phosphoric acid + water"] )
					FRAGMENTS.append( [PREC-MW("PO4H3")-H2O-MW("NH3"), 19, "NL phosphoric acid + water + ammonia"] )
					
					if (self.chains[0][2] > 1):  #if t
						FRAGMENTS.append( [PREC-H2O-H2O, 2, "NL 2xH2O"] )   	
				


			if adduct in NEG_ADDUCTS:
				if adduct == "[M-H]-":  

			
					FRAGMENTS.append([PREC, 199, "precursor"])
					FRAGMENTS.append([MW("PO4H3")-PROTON, 17, "phosphate"])
					FRAGMENTS.append([MW("PO3H")-PROTON, 1000, "phosphite"])

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




x  = LCB_P("LCB_P",[[20,1,1]],  adduct="[M+H]+")
print x.printNist()

