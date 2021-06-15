from lipidRules import *

class LysoCPI(LysoSphingoLipid):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct
			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			
			
			
			if adduct in POS_ADDUCTS:  
				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 55, "precursor"])
					FRAGMENTS.append([PREC-H2O, 18, "NL water"])
					FRAGMENTS.append([PREC-H2O-MW("NH3"), 2, "NL water + ammonia"])			
					FRAGMENTS.append( [PREC-MW("C6H15O10P"),  1000,  "NL inositolphosphate"] )
					FRAGMENTS.append( [PREC-MW("C6H15O10P")+H2O,  120,  "NL deoxyinositolphosphate"] )
					FRAGMENTS.append( [MW("C6H15O10P")-H2O+PROTON,  24,  "deoxyinositol-P"] )														
					
				else:
					FRAGMENTS.append([PREC, 164, "precursor"])
					FRAGMENTS.append([PREC-H2O, 46, "NL water"])
					FRAGMENTS.append([PREC-H2O-MW("NH3"), 7, "NL water + ammonia"])			
					FRAGMENTS.append( [PREC-MW("C6H15O10P")+2*H2O,  49,  "NL doubledeoxyinositolphosphate"] )
					FRAGMENTS.append( [PREC-MW("C6H15O10P")+H2O,  41,  "NL deoxyinositolphosphate"] )									
					FRAGMENTS.append( [MW("C6H12O6")+ADDUCT[adduct],  29,  "inositolphossphate + adduct"] )														
					FRAGMENTS.append( [MW("C6H12O6")-H2O+ADDUCT[adduct],  4,  "deoxyinositol + adduct"] )														
					FRAGMENTS.append( [MW("C6H15O10P")-H2O+ADDUCT[adduct],  1000,  "deoxyinositolphosphate + adduct"] )														
					FRAGMENTS.append( [MW("C6H15O10P")-2*H2O+ADDUCT[adduct],  83,  "deoxyinositol-P + adduct - water"] )														
					
					
			if adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])
			
				FRAGMENTS.append( [MW("PO4H3")-PROTON,  29,  "phosphate"] )
				FRAGMENTS.append( [MW("PO3H")-PROTON,  187,  "phosphite"] )
	   
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON,  25,  "inositolphosphate ion"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O,  563,  "inositolphosphate ion - water"] )	   
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O-H2O,  4,  "inositolphosphate ion - 2x water"] )	   
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H12O6")-PROTON, 21, "NL inositol + adduct"])	   
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-MW("C6H12O6")-PROTON+H2O, 72, "NL deoxyinositol + adduct "])	   	   


					
	

			
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




# x  = LysoCPI("LysoCPI",[[18,1,1]],  adduct="[M-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("LysoCPI", [[18,1,1]], "[M-H]-")