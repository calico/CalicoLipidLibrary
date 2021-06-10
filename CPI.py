from lipidRules import *

class CPI(SphingoLipid):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct
			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]
			FRAGMENTS.append([PREC, 1000, "precursor"])

			
			
			if adduct in NEG_ADDUCTS:  
			
				FRAGMENTS.append( [MW("PO4H3")-PROTON,  1000,  "phosphate"] )
				FRAGMENTS.append( [MW("PO3H")-PROTON,  1000,  "phosphite"] )
	   
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON,  1000,  "phosphoinositol "] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O,  1000,  "phosphinositol - H2O"] )	   
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O-H2O,  1000,  "phosphoinositol - 2xH2O"] )	   
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H12O6")-PROTON, 1000, "Ceramide-P - water"])	   
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H12O6")-PROTON+H2O, 1000, "Ceramide-P"])	   	   
	#LCB fragments:  
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H10O5")-NL(self.chains[1])+MW("CO")-PROTON,1000,"LCB-P - water + CO "])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H10O5")-NL(self.chains[1])-MW("NH")-PROTON,1000,"LCB-P -H3NO "])
	
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H10O5")-NL(self.chains[1])-PROTON,1000,"LCB-P - water "])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H10O5")-NL(self.chains[1])-H2O-PROTON,1000,"LCB-P - 2x water "])
					
	#FA fragments
				FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")-PROTON, 1000, "fatty amide"])
				#FRAGMENTS.append( [NL(self.chains[1])-PROTON, 1000, "fatty acid"])
				if self.chains[1][2] > 0:
					FRAGMENTS.append( [NL(self.chains[1])-MW("CO")-H2O-PROTON, 1000, "fatty acid - CO - H2O"])

			
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




x  = CPI("CPI",[[18,1,1], [24,0,0]],  adduct="[M-H]-")
print x.printNist()

