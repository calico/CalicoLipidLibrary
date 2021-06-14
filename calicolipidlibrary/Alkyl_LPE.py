from lipidRules import *

class Alkyl_LPE(AlkylLysoGPL):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:


				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 540, "precursor"])
					
					FRAGMENTS.append( [PREC-MW("C2H8NO4P"), 9, "NL deoxyethanolamine-P"] )
					FRAGMENTS.append( [PREC-MW("C2H5N"), 7, "NL deoxyethanolamine"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+PROTON, 9, "deoxyethanolamine-P "] )
		
					FRAGMENTS.append( [MW("C2H7NO")+ADDUCT[adduct], 3, "ethanolamine"] )            
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct]-H2O, 31, "glycerol-P - H2O"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct]-H2O, 2, "ethanolamine-P - H2O"])
					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 6, "phosphoric acid"] )
					FRAGMENTS.append( [PREC-H2O, 213, "NL water"] )
					FRAGMENTS.append( [PREC-MW("C3H7O5P"), 415, "NL deoxyglycerol-P"] )
					FRAGMENTS.append( [PREC-MW("C3H7O5P")-H2O, 1000, "NL glycerol-P"] )	

					FRAGMENTS.append( [PREC - MW("C3H6O2"), 70, "NL deoxyglycerol "] )

		
				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 488, "precursor"])
					FRAGMENTS.append( [PREC-MW("C2H8NO4P"), 3, "NL deoxyethanolamine-P"] )
					FRAGMENTS.append( [PREC-MW("C2H6NO3P"), 28, "NL double deoxyethanolamine-P"] )

					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct]-H2O, 270, "glycerol-P - H2O + adduct"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct]-H2O, 46, "double deoxyethanolamine-P + adduct"])
					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 37, "phosphoric acid + Adduct"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct], 3, "deoxyethanolamine-P + adduct"] )
					FRAGMENTS.append( [PREC-MW("C2H5N"), 1000, "NL deoxyethanolamine"] )


			
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])
				FRAGMENTS.append( [MW("H3PO4")-PROTON,  2,  "phosphate"] ) 
				FRAGMENTS.append( [MW("HPO3")-PROTON,  243,  "phosphite"] ) 
				FRAGMENTS.append( [MW("C2H8NO4P")-H2O-PROTON,  10,  "Phosphoethanolamine - H20"] )
				FRAGMENTS.append( [MW("C2H8NO4P")-PROTON,  244,  "Phosphoethanolamine"] )
				FRAGMENTS.append( [MW("C5H12NO5P")-PROTON, 428, "C5H12NO5P"] )
				FRAGMENTS.append( [PREC - MW("C2H7NO"),  51,  " NL ethanolamine"] )

			
					
			
			


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


#calicolipidlibrary.print_spectrum("Alkyl_LPE", [[18,1,0]], "[M+H]+")

# x  = Alkyl_LPE("Alkyl_LPE",[[18,1,0]],  adduct="[M+H]+")
# print x.printNist()

