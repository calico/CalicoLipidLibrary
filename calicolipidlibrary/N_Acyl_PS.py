from lipidRules import *
# needs work
class N_Acyl_PS(NAcylGPL):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]
			FRAGMENTS.append([PREC, 1000, "precursor"])

			if adduct in POS_ADDUCTS:  #NOt done at all
				FRAGMENTS.append( [PREC-(NL(self.chains[2])-MW("O")+MW("NH")+MW("C3H5O6P")), 1000, "NL phosphoserine + Amide"] ) #C3H8NO6P
				FRAGMENTS.append( [NL(self.chains[2])-MW("O")+MW("NH")+MW("C3H2O2")+ADDUCT[adduct], 1000, "Serine+FA+adduct"] )		
					
				if adduct == "[M+H]+":
					#FRAGMENTS.append( [PREC-H2O, 1000, "NL water"] ) #C3H8NO6P
				
					FRAGMENTS.append( [MW("C3H7NO3")+ADDUCT[adduct], 1000,  "Serine"] ) #Ser  
					FRAGMENTS.append( [MW("C3H7NO3")+ADDUCT[adduct]-H2O, 1000,  "Serine - water"] ) #Ser  
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct]-H2O, 1000,  "Glycerol-P - water"] ) 


					for i in [1,0]:
						chain = str(i+1)
# 						FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 1000, "sn" + chain])
# 						FRAGMENTS.append( [PREC - NL(self.chains[i]), 1000, "sn" + chain +"-H2O"])
# 						FRAGMENTS.append( [NL(self.chains[i])-H2O+ADDUCT[adduct], 1000, "sn" + chain +" acylium ion"])
# 						FRAGMENTS.append( [NL(self.chains[i])-H2O-H2O+ADDUCT[adduct], 1000, "sn" + chain +" acylium ion - H2O"])
# 						FRAGMENTS.append( [NL(self.chains[i])-H2O+MW("C3H6O2")+ADDUCT[adduct], 1000, "sn" + chain +" + C3H6O2 (from glycerol)"])
# 						FRAGMENTS.append( [NL(self.chains[i])-H2O+MW("C3H5NO2")+ADDUCT[adduct], 1000, "sn" + chain +" + C3H5NO2"])




				if adduct == "[M+Na]+" or adduct == "[M+NH4]+" or adduct == "[M+K]+":
					FRAGMENTS.append( [PREC-(NL(self.chains[2])-MW("O")+MW("NH")+MW("C3H5O6P"))-ADDUCT[adduct]+PROTON, 1000, "NL phosphoserine + Amide - adduct"] ) #C3H8NO6P				
					FRAGMENTS.append( [NL(self.chains[2])-MW("O")+MW("NH")+MW("C3H5O6P")+ADDUCT[adduct], 1000, "PhosphoSerine+FA"] )		
					FRAGMENTS.append( [NL(self.chains[2])-MW("O")+MW("NH")+MW("C3H2O2")+PROTON, 1000, "Serine+FA + proton"] )		
					FRAGMENTS.append( [PREC - (NL(self.chains[2])-MW("O")+MW("NH")+MW("C3H2O2")), 1000, "NL Serine+FA"] )		

					#FRAGMENTS.append( [PREC-MW("C3H8NO6P")+PROTON-ADDUCT[adduct], 1000,  "precursor - P-Ser - adduct"] ) 
					#FRAGMENTS.append( [MW("C3H8NO6P")+ADDUCT[adduct], 1000,  " P-Ser+adduct"] )
					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 1000,  " H3PO4+adduct"] )

	
	
	

		
			elif adduct in NEG_ADDUCTS:  #NEG is done, at least for M-H
		
				FRAGMENTS.append( [MW("H3PO4")-PROTON,  1000,  "H2PO4-"] ) 
				FRAGMENTS.append( [MW("HPO3")-PROTON,  1000,  "PO3-"] ) 
	
				FRAGMENTS.append( [MW("C3H7O5P")-PROTON, 1000, "C3H6O5P-"] )
				FRAGMENTS.append( [MW("C3H9O6P")-PROTON, 1000, "C3H8O6P-"] )
				FRAGMENTS.append( [PREC - (NL(self.chains[2])-MW("O")+MW("NH")+MW("C3H4O2")), 1000, "NL Serine+FA"] )		

				for i in [0,1]:
					chain = str(i+1)
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn" + chain +" RCO2-"])
					FRAGMENTS.append( [NL(self.chains[i])+MW("C3H9O6P")-H2O-PROTON, 1000, "sn" + chain +" Glycerol-P"])
					FRAGMENTS.append( [NL(self.chains[i])+MW("C3H9O6P")-2*H2O-PROTON, 1000, "sn" + chain +" Glycerol-P-H2O"])
					#FRAGMENTS.append( [PREC - MW("C3H5NO2")-NL(self.chains[i])+H2O, 1000, "NL (Serine + sn" + chain +")"] )
					#FRAGMENTS.append( [PREC - MW("C3H5NO2")-NL(self.chains[i]), 1000, "NL Serine + sn" + chain +" + H2O)"] )

	

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



#
# x  = N_Acyl_PS("N_Acyl_PS",[[16,0,0], [18,1,0], [18,2,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("N_Acyl_PS",[[16,0,0], [18,1,0], [18,2,0]], "[M+H]+")