from lipidRules import *

class CPE(SphingoLipid):





		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
				


			
				if adduct == "[M+H]+":  
					FRAGMENTS.append([PREC, 6, "precursor"])
					FRAGMENTS.append([PREC-H2O, 90, "NL water"])

					FRAGMENTS.append( [PREC - MW("C2H7NO"), 47, "NL ethanolamine"] )       
					FRAGMENTS.append( [MW("C2H10NO5P")-H2O+ADDUCT[adduct], 6, "deoxyethanolamine phosphate"] )       
					FRAGMENTS.append( [PREC - MW("C2H10NO5P") , 680, "NL phosphoethanolamine"] )
					FRAGMENTS.append( [PREC - MW("C2H10NO5P")  + H2O, 1000, "NL deoxyphosphoethanolamine"] )
					
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("C2H3N")+PROTON, 134, "Fatty amide + C2"])

					FRAGMENTS.append( [PREC - NL(self.chains[1]) +H2O - MW("C2H10NO5P") - ADDUCT[adduct] + PROTON, 455, "NL phosphoethanolamine + ketene + ADDUCT"] ) 
	
	
	
				if adduct == "[M+Na]+" or adduct == "[M+K]+" or adduct == "[M+NH4]+":
					FRAGMENTS.append([PREC, 787, "precursor"])

					FRAGMENTS.append( [PREC - MW("C2H7NO"), 19, "NL ethanolamine"] )       
					FRAGMENTS.append( [PREC - MW("C2H7NO")+H2O, 318, "NL deoxyethanolamine"] )       
					FRAGMENTS.append( [MW("C2H10NO5P")-H2O+ADDUCT[adduct], 235, "deoxyethanolamine phosphate + Adduct"] )       
					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 189, "phosphoric acid + Adduct"] )       
					FRAGMENTS.append( [PREC - MW("C2H10NO5P") , 25, "NL phosphoethanolamine"] )
					FRAGMENTS.append( [PREC - MW("C2H10NO5P")  + H2O, 955, "NL deoxyhosphoethanolamine"] )

					
					FRAGMENTS.append( [PREC - MW("C2H10NO5P") - (ADDUCT[adduct]-PROTON), 57, "NL phosphoethanolamine + adduct"] )
					FRAGMENTS.append( [PREC - MW("C2H10NO5P") - (ADDUCT[adduct]-PROTON) + H2O, 25, "NL Deoxyphosphoethanolamine + adduct"] )
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("C2H3N")+PROTON, 16, "Fatty amide + C2"])

					FRAGMENTS.append( [PREC - NL(self.chains[1]) +H2O - MW("C2H10NO5P") - ADDUCT[adduct] + PROTON, 23, "NL Phosphoethanolamine + ketene + ADDUCT"] ) 
	
			
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])

				FRAGMENTS.append( [MW("C3H9O6P")-H2O-PROTON, 1000,  "glycerol-P - H2O"] ) 
				FRAGMENTS.append( [MW("H3PO4")-PROTON,  2,  "phosphate"] ) 

				FRAGMENTS.append( [MW("C2H8NO4P")-PROTON,  385,  "phosphoethanolamine"] )

				FRAGMENTS.append( [PREC - MW("C2H7NO"),  2,  "NL ethanolamine"] )
				FRAGMENTS.append( [PREC - MW("C2H7NO")+H2O,  2,  "NL deoxyethanolamine"] )

				i=1
				chain = str(i+1)
				FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 7, "NL ketene + 2x water"])
		
			
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



# x  = CPE("CPE",[[17,1,1], [12,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# x  = CPE("CPE",[[18,1,1], [24,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# x  = CPE("CPE",[[18,1,1], [24,1,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("CPE", [[18,1,1], [24,1,0]], "[M+H]+")