from lipidRules import *

class Alkyl_LPC(AlkylLysoGPL):


		neg_adduct_set = GPL.neg_adduct_set + ["[M-CH3]-"]



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:


				if adduct == "[M+H]+":  
					FRAGMENTS.append([PREC, 1000, "precursor"])

					FRAGMENTS.append( [MW("C4H8N")+PROTON, 5, "C4H9N+"] )
					FRAGMENTS.append( [MW("C5H13NO")-H2O+PROTON, 127, "Dehydro Choline"] )
					FRAGMENTS.append( [MW("C3H9N")+PROTON, 15, "Trimethylamine"] )
					FRAGMENTS.append( [MW("C8H20NO6P")-H2O+PROTON, 108, "Glycero-phosphocholine - water"] )     
					FRAGMENTS.append( [MW("C8H20NO6P")-MW("C3H9N")-H2O+PROTON, 101, "Glycero-phosphocholine - water - trimethylamine"] )        	
					FRAGMENTS.append( [MW("C5H14NO4P")+ADDUCT[adduct], 180, "Phospocholine"] )      	
					FRAGMENTS.append( [MW("C2H5O4P")+PROTON, 80, "C2H5O4P + PROTON"] )
					FRAGMENTS.append( [MW("C5H13NO")+ADDUCT[adduct], 417, "Choline"] )
					FRAGMENTS.append( [MW("H3O4P")+ADDUCT[adduct], 1, "phosphoric acid"] )

		
				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 744, "precursor"])

					FRAGMENTS.append( [MW("C4H8N")+PROTON, 5, "C4H8N + PROTON"] )
					FRAGMENTS.append( [MW("C5H13NO")-H2O+PROTON, 113, "Dehydro Choline"] )
					FRAGMENTS.append( [MW("C3H9N")+PROTON, 5, "Trimethylamine"] )
					FRAGMENTS.append( [MW("C8H20NO6P")-MW("C3H9N")-H2O+ADDUCT[adduct], 2, "Glycero-phosphocholine - water - trimethylamine + adduct"] )        	
					FRAGMENTS.append( [PREC-MW("C5H13NO"), 11, "NL choline"] )                	
					FRAGMENTS.append( [MW("C2H5O4P")+ADDUCT[adduct], 1000, "C2H5O4P + Adduct"] )
					FRAGMENTS.append( [MW("C5H13NO")+PROTON, 712, "Choline"] )
					FRAGMENTS.append( [PREC - MW("C3H9N"), 1000, "NL Trimethylamine"] )      
	 


			
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 100, "precursor"])
				FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON-MW("CH3"), 1000, "NL CH3 + Adduct"] )
				FRAGMENTS.append( [MW("H3PO4")-PROTON,  50,  "phosphate"] )
				FRAGMENTS.append( [MW("HPO3")-PROTON,  200,  "phosphite"] )

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


# calicolipidlibrary.print_spectrum("Alkyl_LPC",[[18,1,0]], "[M+H]+")

# x  = Alkyl_LPC("Alkyl_LPC",[[18,1,0]],  adduct="[M+H]+")
# print x.printNist()
#
# x  = Alkyl_LPC("Alkyl_LPC",[[18,1,0]],  adduct="[M+Li]+")
# print x.printNist()


