from lipidRules import *

class DGDG(GPL):

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:


			
				if adduct == "[M+Na]+" or adduct == "[M+K]+":
					FRAGMENTS.append([PREC, 1000, "precursor"])
	
					FRAGMENTS.append( [PREC-MW("C6H10O5"), 1000, "NL deoxyhexose"] )

					FRAGMENTS.append( [MW("C9H16O6")+MW("C6H10O5")+ADDUCT[adduct], 1000, "dihexosylglycerol - 2 water +Adduct"] )
					FRAGMENTS.append( [MW("C12H22O11")-H2O+ADDUCT[adduct], 1000, "dihexose - water + Adduct"] )



				for i in range(0,len(self.chains)):
						chain = str(i+1)
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O+MW("C3H6O2"), 1000, "sn"+chain+" acylium ion + glycerol- water"] )
						FRAGMENTS.append( [PREC - NL(self.chains[i]), 1000, "NL sn"+chain] )
				
				if adduct == "[M+NH4]+":  
					FRAGMENTS.append([PREC, 1000, "precursor"])
	
					
					FRAGMENTS.append( [PREC - MW("C12H22O11")+H2O-ADDUCT[adduct]+PROTON, 1000, "NL deoxydihexose + adduct"] )
					FRAGMENTS.append( [PREC - MW("C12H22O11")-ADDUCT[adduct]+PROTON, 1000, "NL dihexose +adduct "] )


					for i in range(0,len(self.chains)):
						chain = str(i+1)
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O+MW("C3H6O2"), 1000, "sn"+chain+" acylium ion + glycerol - water"] )



 
			if adduct in NEG_ADDUCTS:
				FRAGMENTS.append([PREC, 1000, "precursor"])
				for i in range(0,len(self.chains)):
					chain = str(i+1)
					FRAGMENTS.append([NL(self.chains[i])-PROTON, 1000, "sn" + chain])	
					FRAGMENTS.append([PREC - NL(self.chains[i]), 1000, "NL sn" + chain])	
					FRAGMENTS.append([PREC - NL(self.chains[i])+H2O, 1000, "NL sn" + chain + " ketene"])	

				FRAGMENTS.append([PREC - NL(self.chains[0])-NL(self.chains[1])+2*H2O, 1000, "diglucose + glycerol"])	
				FRAGMENTS.append([PREC - NL(self.chains[0])-NL(self.chains[1])+H2O, 1000, "diglucose + deoxyglycerol"])	
				FRAGMENTS.append([PREC - NL(self.chains[0])-NL(self.chains[1]), 1000, "diglucose +deoxyglycerol - water"])	
				


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




# x  = DGDG("DGDG",[[18,3,0], [16,3,0]],  adduct="[M+NH4]+")
# print x.printNist()
#
#
#
# x  = DGDG("DGDG",[[18,3,0], [16,3,0]],  adduct="[M+Na]+")
# print x.printNist()
#
#
#
# x  = DGDG("DGDG",[[18,3,0], [16,3,0]],  adduct="[M-H]-")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("DGDG", [[18,3,0], [16,3,0]], "[M+NH4]+")