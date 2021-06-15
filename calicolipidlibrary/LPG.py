from lipidRules import *

class LPG(LysoGPL):

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]


			if adduct in POS_ADDUCTS:  


				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 24, "precursor"])				
					FRAGMENTS.append( [PREC - MW("C3H9O6P"), 1000, "NL glycerol-P"])

					FRAGMENTS.append([PREC-2*H2O, 64, "precursor - 2x water"])
					FRAGMENTS.append([PREC-MW("C3H8O3"), 5, "precursor - glycerol"])
										
					FRAGMENTS.append([MW("C3H9O6P")+ADDUCT[adduct], 18, "glycerol-P + adduct"])
					FRAGMENTS.append([MW("C3H9O6P")-H2O+ADDUCT[adduct], 72, "deoxyglycerol-P + adduct"])
					
					FRAGMENTS.append( [NL(self.chains[0])-H2O+PROTON, 26, "RCO+ sn1"])
					FRAGMENTS.append( [PREC - NL(self.chains[0]), 4, "NL sn1"])
					FRAGMENTS.append( [MW("H3O4P")+ADDUCT[adduct], 25, "phosphoric acid + adduct"])

					FRAGMENTS.append( [NL(self.chains[0])-H2O+PROTON, 26, "sn1 acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[0])+H2O, 4, "NL ketene"])


				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 976, "precursor"])
					FRAGMENTS.append([PREC-H2O, 49, "NL  water"])
					FRAGMENTS.append([PREC-MW("C3H6O2"), 49, "NL deoxyglycerol"])
					FRAGMENTS.append([PREC-MW("C3H8O3"), 49, "NL glycerol"])
					FRAGMENTS.append([PREC-MW("C3H9O6P"), 49, "NL glycerol-P"])
										
					FRAGMENTS.append([MW("C3H9O6P")+ADDUCT[adduct]+H2O, 8, "glycerol-P + water + adduct"])								
					FRAGMENTS.append([MW("C3H9O6P")+ADDUCT[adduct], 281, "glycerol-P + adduct"])
					FRAGMENTS.append([MW("C6H15O8P")+ADDUCT[adduct], 4, "glycerol-P-glycerol + adduct"])

					FRAGMENTS.append([MW("C3H9O6P")-H2O+ADDUCT[adduct], 776, "deoxyglycerol-P + adduct"])

					FRAGMENTS.append([MW("C3H8O3")+ADDUCT[adduct], 14, "glycerol + adduct"])
					
					FRAGMENTS.append( [NL(self.chains[0])-H2O+PROTON, 8, "sn1 acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[0]), 251, "NL sn1"])
					FRAGMENTS.append( [MW("H3O4P")+ADDUCT[adduct], 9, "phosphoric acid + adduct"])

			
			
			elif adduct in NEG_ADDUCTS:

				if self.chains[0][0] != 0:
					i = 0
				else:
					i = 1


				FRAGMENTS.append([PREC, 385, "precursor"])

				FRAGMENTS.append( [MW("H3PO4")-PROTON,  13,  "phosphate"] ) 
				FRAGMENTS.append( [MW("HPO3")-PROTON,  1,  "phosphite"] ) 
				FRAGMENTS.append( [MW("C3H7O5P")-PROTON, 153, "glycerol-P - water"] )
				FRAGMENTS.append( [MW("C3H9O6P")-PROTON, 10, "glycerol-P"] )

				FRAGMENTS.append( [PREC - MW("C3H7O3")-PROTON, 6, "NL glycerol"] )
				FRAGMENTS.append( [PREC - MW("C3H7O3")+H2O-PROTON, 2, "NL deoxyglycerol"] )

				FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn1"])
				FRAGMENTS.append( [PREC-NL(self.chains[i]), 47, "NL sn1"])
				FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O, 12, "NL sn1 ketene "])

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




# x  = LPG("LPG",[[18,1,0],[0,0,0]],  adduct="[M+Na]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("LPG", [[18,1,0],[0,0,0]], "[M+Na]+")
