from lipidRules import *

class Alkyl_PE(alkylGPL):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
			
			
				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 23, "precursor"])
					FRAGMENTS.append( [PREC-H2O, 12, "NL water"] )
					FRAGMENTS.append( [PREC-MW("C2H8NO4P"), 39, "NL ethanolamine-P"] )
			
					FRAGMENTS.append( [PREC - NL(self.chains[1])+H2O, 9, "NL sn2 ketene"])
					FRAGMENTS.append( [PREC - NL(self.chains[1]), 4, "NL sn2"])
					FRAGMENTS.append( [PREC-NL(self.chains[1])+H2O - MW("C3H6O2"), 345, "NL sn2 + glycerol"])            		    		
			
					FRAGMENTS.append( [PREC-MW("C2H8NO4P")-ADDUCT[adduct]-NL(self.chains[1])+PROTON, 31, "NL sn2 + ethanolamine-P"] ) 
					FRAGMENTS.append( [PREC-MW("C3H7O5P")-ADDUCT[adduct]-NL(self.chains[1])+PROTON, 109, "NL sn2 + glycerol-phosphate"] ) 

					FRAGMENTS.append( [NL(self.chains[1])-H2O+ MW("C3H6O2")+ADDUCT[adduct], 1000, "sn2 + deoxyglycerol"])            		    				

		
				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 10000, "precursor"])

					FRAGMENTS.append( [PREC-MW("C2H8NO4P"), 304, "NL ethanolamine-P"] )
					FRAGMENTS.append( [PREC-MW("C2H8NO4P")+H2O, 54, "NL deoxyethanolamine-P"] )
				
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct]-H2O, 27, "glycerol-P - H2O"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct]-H2O, 43, "ethanolamine-P - H2O"])
					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 168, "phosphoric acid + adduct"] )				
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct], 209, "ethanolamine-P + adduct"] )
					FRAGMENTS.append( [PREC-MW("C2H5N"), 559, "NL deoxyethanolamine"] )
					FRAGMENTS.append( [PREC-MW("C2H8NO4P")+PROTON-ADDUCT[adduct], 26, "NL ethanolamine phosphate + adduct"] )

					FRAGMENTS.append( [NL(self.chains[0])-H2O-H2O+PROTON, 1, "sn1 ether-linked hydrocarbon"])  # is this from sn1 or sn2?  

					
					FRAGMENTS.append( [PREC-NL(self.chains[1])+H2O - MW("C3H6O2"), 21, "NL sn2 + glycerol"])            		    		
					FRAGMENTS.append( [PREC-MW("C2H5N")-NL(self.chains[1]), 14, "NL sn2 + deoxyethanolamine"] ) 
					FRAGMENTS.append( [PREC-MW("C2H8NO4P")-ADDUCT[adduct]-NL(self.chains[1])+PROTON, 78, "NL sn2 + ethanolamine-phosphate"] ) 
					FRAGMENTS.append( [PREC-NL(self.chains[0])+H2O-MW("H4")-MW("C2H5N"), 115, "NL sn1 alcohol + deoxyethanolamine"])            		    		
		

			
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])

				FRAGMENTS.append( [MW("C2H8NO4P")-PROTON,  70,  "Phosphoethanolamine"] )
				FRAGMENTS.append( [MW("C5H12NO5P")-PROTON, 83, "phosphoethanolamine+glycerol"] )

				FRAGMENTS.append( [NL(self.chains[1])-PROTON, 1000, "sn2"])
				FRAGMENTS.append( [PREC - NL(self.chains[1]) +H2O - ADDUCT[adduct] - PROTON, 369, "NL sn2 ketene"])
				FRAGMENTS.append( [PREC - NL(self.chains[1]) - ADDUCT[adduct] - PROTON, 124, "NL sn2 - water"])
				FRAGMENTS.append( [NL(self.chains[0])-H2O+MW("H4")-PROTON, 1, "sn1 ether-linked alcohol"]) 
					
				if self.chains[1][1] > 3:   #if have PUFA SN2
					FRAGMENTS.append( [NL(self.chains[1])-MW("CO2")-PROTON, 409, "sn2 - CO2"])


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




x  = Alkyl_PE("Alkyl_PE",[[18,1,0], [22,6,0]],  adduct="[M-H]-")
print x.printNist()

