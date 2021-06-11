from lipidRules import *

class MIP2C(SphingoLipid):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

	
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]
			FRAGMENTS.append([PREC, 1000, "precursor"])

		
		
			if adduct in NEG_ADDUCTS:  
		
				FRAGMENTS.append( [MW("C18H33O19P")-PROTON,  1000,  "MIP2 - PO3 "] )	   
				FRAGMENTS.append( [MW("C12H23O14P")-PROTON,  1000,  "MIP "] )
				FRAGMENTS.append( [MW("C12H23O14P")-PROTON-H2O,  1000,  "MIP - H2O"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON,  1000,  "phosphoinositol ion"] )	   
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O,  1000,  "phosphoinositol ion-H2O"] )	   
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O-H2O,  1000,  "phosphoinositol ion-2xH2O"] )	
			   
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C18H31O18P")-PROTON-H2O, 1000, "Ceramide-P - H2O "])	   
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C18H31O18P")-PROTON, 1000, "Ceramide-P"])	 
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C18H31O18P")-MW("PO3H")-PROTON, 1000, "Ceramide"])	   	   
				   
				#LCB fragments:  
			
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-NL(self.chains[1])-PROTON,1000,"LCB-P - water "])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-NL(self.chains[1])-H2O-PROTON,1000,"LCB-P -2x water "])

				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-NL(self.chains[1])+MW("CO")-PROTON,1000,"LCB-P - water  + CO "])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-NL(self.chains[1])-MW("NH")-PROTON,1000,"LCB-P - H3NO "])

			
#				FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")-PROTON, 1000, "sn2 fatty amide"])
				FRAGMENTS.append( [NL(self.chains[1])-PROTON, 1000, "fatty acid"])  #seems wrong?
				if self.chains[1][2] > 0:
					FRAGMENTS.append( [NL(self.chains[1])-MW("CO")-H2O-PROTON, 1000, "fatty acid - CO - H2O"])

		
			return(FRAGMENTS)


		@classmethod
		def generateLibrary(cls,target=None, mode="pos"):
			if target: handle = open(target, 'a+')
			if mode == "pos":  
				adduct_set = self.pos_adduct_set
			elif mode =="neg":
				adduct_set = self.neg_adduct_set
			for c1 in self.chain1_ranges:
				for c2 in self.chain2_ranges:
					for adduct in adduct_set:
						x  = MIP2C("MIP2C",[c1[0],c2[0]], [c1[1],c2[1]], [c1[2],c2[2]],  adduct=adduct)						
						content = x.printNist()
						if target: handle.write(content)
						else: sys.stdout.write(content)
			if target: handle.close()



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




x  = MIP2C("MIP2C",[[18,1,1], [24,0,0]],  adduct="[M+H]+")
print x.printNist()

