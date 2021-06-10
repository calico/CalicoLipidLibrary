from lipidRules import *

class LPS(LysoGPL):




		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
		
			
				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 19, "precursor"])
				
					FRAGMENTS.append( [PREC-H2O, 47, "NL water"] ) #C3H8NO6P
					FRAGMENTS.append( [PREC-MW("C3H8NO6P"), 1000, "NL phosphoserine"] ) #C3H8NO6P
					FRAGMENTS.append( [MW("C3H8NO6P")+ADDUCT[adduct], 6,  "phosphoserine"] ) #Ser       	
					FRAGMENTS.append( [MW("C3H7NO3")+ADDUCT[adduct], 178,  "Serine"] ) #Ser  
					FRAGMENTS.append( [MW("C3H7NO3")+ADDUCT[adduct]-H2O, 215,  "Serine - water"] ) #Ser  
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct], 43,  "Glycerol-P"] ) #Ser  C3H9O6P
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct]-H2O, 151,  "Glycerol-P - water"] ) 

					FRAGMENTS.append( [PREC-MW("C3H7NO3"), 309, "NL serine"] ) #C3H8NO6P					
	
					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 141,  "phosphoric acid + Adduct"] ) #Ser       				
					FRAGMENTS.append( [PREC - NL(self.chains[0])+H2O, 18, "NL sn1 ketene"])
					FRAGMENTS.append( [PREC - NL(self.chains[0]), 10, "NL sn1"])
					FRAGMENTS.append( [NL(self.chains[0])-H2O+ADDUCT[adduct], 99, "sn1 acylium ion"])
					FRAGMENTS.append( [NL(self.chains[0])-H2O-H2O+ADDUCT[adduct], 25, "sn1 acylium ion - H2O"])
					FRAGMENTS.append( [NL(self.chains[0])-H2O+MW("C3H5NO2")+ADDUCT[adduct], 85, "sn1 + serine"])


				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:

					FRAGMENTS.append([PREC, 1000, "precursor"])

					FRAGMENTS.append( [PREC-H2O, 41, "NL water"] ) #C3H8NO6P
					FRAGMENTS.append( [PREC-MW("C3H8NO6P"), 16, "NL phosphoserine"] ) #C3H8NO6P
					FRAGMENTS.append( [PREC-MW("C3H8NO6P")-ADDUCT[adduct]+PROTON, 439, "NL phosphoserine + adduct"] ) #C3H8NO6P
					
					FRAGMENTS.append( [PREC-MW("C3H8NO6P")+H2O, 49, "NL deoxyphosphoserine"] ) #C3H8NO6P

					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 141,  "Phosphoric acid + Adduct"] ) #Ser       				
					FRAGMENTS.append( [MW("C3H8NO6P")+ADDUCT[adduct], 207,  "phosphoserine + Adduct"] ) #Ser       	
					FRAGMENTS.append( [MW("C3H7NO3")+ADDUCT[adduct], 6,  "Serine + Adduct"] ) #Ser  
					FRAGMENTS.append( [MW("C3H7NO3")+ADDUCT[adduct]-H2O, 25,  "Serine - water + Adduct"] ) #Ser  
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct], 79,  "Glycerol-P + Adduct"] ) #Ser  C3H9O6P
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct]-H2O, 45,  "Glycerol-P - water + Adduct"] ) 

					FRAGMENTS.append( [PREC-MW("C3H7NO3"), 309, "NL serine"] ) #C3H8NO6P					
					FRAGMENTS.append( [PREC-MW("C3H5NO2"), 887, "NL deoxyserine"] ) #C3H8NO6P					


					FRAGMENTS.append( [PREC - NL(self.chains[0])+H2O, 21, "NL sn1 ketene"])
					FRAGMENTS.append( [PREC - NL(self.chains[0]), 1, "NL sn1"])

			
			
			
			elif adduct in NEG_ADDUCTS:

				if self.chains[0][0] != 0:
					i = 0
				else:
					i = 1


				FRAGMENTS.append([PREC, 121, "precursor"])
			
				FRAGMENTS.append( [MW("H3PO4")-PROTON, 11,  "phosphate"] ) 
				FRAGMENTS.append( [MW("HPO3")-PROTON,  142,  "phosphite"] ) 
				FRAGMENTS.append( [PREC - MW("C3H5NO2"), 443, "NL serine"] )
				FRAGMENTS.append( [PREC - MW("C3H5NO2")-H2O, 3, "NL serine + water"] )
		
				FRAGMENTS.append( [MW("C3H7O5P")-PROTON, 1000, "phosphoserine - water"] )
				FRAGMENTS.append( [MW("C3H9O6P")-PROTON, 37, "phosphoserine"] )
		
				chain = str(i+1)
				FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn" + chain])


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




x  = LPS("LPS",[[18,1,0],[0,0,0]],  adduct="[M+H]+")
print x.printNist()

