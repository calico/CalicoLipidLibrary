from lipidRules import *

class LysoCPE(LysoSphingoLipid):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
				


			
				if adduct == "[M+H]+":  
					FRAGMENTS.append([PREC, 16, "precursor"])
					FRAGMENTS.append([PREC-H2O, 15, "NL water"])
					FRAGMENTS.append([PREC-H2O-MW("NH3"), 22, "NL water + ammonia"])
					FRAGMENTS.append( [PREC - MW("C2H5N"), 17, "NL deoxyethanolamine"])

					FRAGMENTS.append( [PREC - MW("C2H7NO"), 21, "NL ethanolamine"] )       
					FRAGMENTS.append( [MW("C2H10NO5P")-H2O+ADDUCT[adduct], 72, "deoxyethanolamine phosphate"] )       
					FRAGMENTS.append( [PREC - MW("C2H10NO5P") , 1000, "NL Phosphoethanolamine"] )
					FRAGMENTS.append( [PREC - MW("C2H10NO5P")  + H2O, 92, "NL Deoxyphosphoethanolamine"] )
					
					FRAGMENTS.append( [ MW("C2H7NO)")+PROTON, 17, "ethanolamine"])

					#additional peak at 82.06, but I cannot figure out.  
	
	
	
				if adduct == "[M+Na]+" or adduct == "[M+K]+":  
					FRAGMENTS.append([PREC, 787, "precursor"])
					FRAGMENTS.append([PREC-H2O, 90, "NL water"])
					FRAGMENTS.append([PREC-H2O-MW("NH3"), 4, "NL water + ammonia"])

					FRAGMENTS.append( [PREC - MW("C2H7NO"), 35, "NL ethanolamine"] )       
					FRAGMENTS.append( [PREC - MW("C2H7NO")+H2O, 100, "NL deoxyethanolamine"] )       

					  
					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 600, "phosphoric acid + adduct"] )       
					FRAGMENTS.append( [PREC - MW("C2H10NO5P") , 7, "NL phosphoethanolamine"] )
					FRAGMENTS.append( [PREC - MW("C2H10NO5P")  + H2O, 955, "NL deoxyhosphoethanolamine"] )
					FRAGMENTS.append( [MW("C2H10NO5P")-2*H2O+ADDUCT[adduct], 24, "doubledeoxyethanolamine phosphate"] )       
					FRAGMENTS.append( [MW("C2H10NO5P")-H2O+ADDUCT[adduct], 1000, "deoxyethanolamine phosphate + adduct"] )       

					
					FRAGMENTS.append( [PREC - MW("C2H10NO5P") - (ADDUCT[adduct]-PROTON), 37, "NL Phosphoethanolamine + adduct"] )
	
			
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])

				FRAGMENTS.append( [MW("H3PO4")-PROTON,  11,  "phosphate"] ) 
				FRAGMENTS.append( [MW("HPO3")-PROTON,  550,  "phosphite"] ) 

				FRAGMENTS.append( [MW("C2H8NO4P")-PROTON,  578,  "Phosphoethanolamine"] )

				FRAGMENTS.append( [PREC - MW("C2H7NO"),  9,  "NL ethanolamine"] )
				FRAGMENTS.append( [PREC - MW("C2H7NO")+H2O,  12,  "NL deoxyethanolamine"] )
		
				#additional peak at 181 I cannot figure out
			
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



x  = LysoCPE("LysoCPE",[[18,1,1]],  adduct="[M-H]-")
print x.printNist()

