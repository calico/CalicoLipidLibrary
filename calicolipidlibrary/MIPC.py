from lipidRules import *

class MIPC(SphingoLipid):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

#			print self.hydroxyls[0]  #DEBUG line...seems to be one too high, given how I define hydroxyls on sbases
#			print self.hydroxyls[1]
		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]
			FRAGMENTS.append([PREC, 1000, "precursor"])


			if adduct in NEG_ADDUCTS:  
		
				FRAGMENTS.append( [MW("PO4H3")-PROTON,  88,  "phosphate"] )
				FRAGMENTS.append( [MW("PO3H")-PROTON,  607,  "phosphite"] )
   
				FRAGMENTS.append( [MW("C12H23O14P")-PROTON,  330,  "MIP "] )
				FRAGMENTS.append( [MW("C12H23O14P")-PROTON-H2O,  95,  "MIP - H2O"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON,  15,  "inositolphosphate ion"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O,  25,  "inositolphosphate ion-H2O"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O-H2O,  8,  "inositolphosphate ion-2xH2O"] )
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-PROTON-H2O, 26, "Ceramide-P - H2O "])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-PROTON, 10, "Ceramide-P"])
				#LCB fragments:  
			
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-NL(self.chains[1])-PROTON,1,"LCB-P - water "])
				#FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-NL(self.chains[1])-H2O-PROTON,1000,"LCB-P - 2x water "])

				#FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-NL(self.chains[1])+MW("CO")-PROTON,1000,"LCB-P - water + CO "])
				#FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C12H20O10")-NL(self.chains[1])-MW("NH")-PROTON,1000,"LCB-P - H3NO "])

			
#				FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")-PROTON, 1000, "sn2 fatty amide"])
				#FRAGMENTS.append( [NL(self.chains[1])-PROTON, 1000, "fatty amide"])

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




x  = MIPC("MIPC",[[18,0,2], [26,0,0]],  adduct="[M-H]-")
print x.printNist()

x  = MIPC("MIPC",[[18,0,1], [26,0,1]],  adduct="[M-H]-")
print x.printNist()
