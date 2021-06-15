from lipidRules import *

class LCB(LysoSphingoLipid):

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
				FRAGMENTS.append([PREC, 180, "precursor"])

				if adduct == "[M+H]+":  
					FRAGMENTS.append( [PREC-H2O, 1000, "NL water"] )
					FRAGMENTS.append( [PREC-H2O-H2O-MW("NH3"), 4, "NL 2x water + ammonia"] )   						        	
					
					if (self.chains[0][2] > 0):  #if d or t
						FRAGMENTS.append( [PREC-H2O-H2O, 136, "NL 2 x water"] )   	
					if (self.chains[0][1] == 0): #if has no double bond
						FRAGMENTS.append( [PREC-MW("CH4O2N"), 47, "NL CH4O2N"] )   	        	
					FRAGMENTS.append( [PREC-MW("CH4O2"), 214, "NL CH4O2"] )   	        	
					FRAGMENTS.append( [MW("C2H5NO")+PROTON, 662, "C2H6ON"] )   	        	





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




# x  = LCB("LCB",[[18,1,1]],  adduct="[M+H]+")
# print x.printNist()
#
# x  = LCB("LCB",[[18,0,1]],  adduct="[M+H]+")
# print x.printNist()
#
# x  = LCB("LCB",[[20,0,1]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("LCB", [[18,1,1]], "[M+H]+")