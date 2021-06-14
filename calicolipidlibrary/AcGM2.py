from lipidRules import *

class AcGM2(SphingoLipid):

		pos_adduct_set = ["[M+H]+"]
		neg_adduct_set = ["[M-H]-"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 1, "precursor"])

					FRAGMENTS.append( [PREC-H2O, 38, "NL water"] )

					FRAGMENTS.append( [PREC-MW("C11H19NO9"), 62, "NL N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C8H15NO6"), 23, "NL N-Ac-Glucosamine"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")+H2O, 76, "NL N-Ac-Neuraminic Acid - N-Ac-Glucosamine "] )
										
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H10O5")+H2O, 76, "NL N-Ac-Neuraminic Acid - N-Ac-Glucosamine - glucose"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")+H2O, 255, "Ceramide"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5"), 54, "Ceramide - water"] )

					FRAGMENTS.append( [PREC - MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")-NL(self.chains[1])+H2O*2, 71, "sphingoid base"])
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")-NL(self.chains[1])+H2O, 616, "sphingoid base - water"])
					FRAGMENTS.append( [PREC - MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")-NL(self.chains[1])+H2O-MW("C"), 56, "sphingoid base - CH2O"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O+MW("NH3")+PROTON, 16, "sn2 fatty amide"])

					FRAGMENTS.append([MW("C11H19NO9") + PROTON, 9, "N-Ac-Neuraminic Acid"])
					FRAGMENTS.append([MW("C11H19NO9")-H2O + PROTON, 102, "deoxy-N-Ac-Neuraminic Acid"])

					FRAGMENTS.append([MW("C8H15NO6") + PROTON, 19, "N-Ac-glucosamine "])
					FRAGMENTS.append([MW("C8H15NO6") - H2O + PROTON, 1000, "N-Ac-glucosamine - water"])
					FRAGMENTS.append([MW("C8H15NO6") - 2 * H2O + PROTON, 461, "N-Ac-glucosamine - 2*water"])




			# if adduct == "[M+Na]+" or adduct == "[M+K]+":
			# 		FRAGMENTS.append( [PREC - MW("CO2"), 1000, "Precursor - CO2"] )
			# 		#FRAGMENTS.append( [MW("C11H19NO10")+ADDUCT[adduct], 1000, "pre - N-Ac-Neuraminic Acid"] )
			# 		FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 1000, "pre - N-Ac-Neuraminic Acid"] )
			# 		FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")+2*H2O, 1000, "pre - N-Ac-Neuraminic Acid - N-Ac-Glucosamine"] )
			# 		FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H10O5")+2*H2O, 1000, "pre - N-Ac-Neuraminic Acid - N-Ac-Glucosamine - glucose"] )
			# 		FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")+2*H2O, 1000, "pre - N-Ac-Neuraminic Acid - N-Ac-Glucosamine - glucose - galactose"] )
			# 		FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")+H2O, 1000, "pre - N-Ac-Neuraminic Acid - N-Ac-Glucosamine - glucose - galactose - water"] )
			# 		FRAGMENTS.append( [PREC - MW("C11H19NO9")-2*MW("C6H10O5")-NL(self.chains[1])+H2O, 1000, "pre- (head + sn2)"])
			# 		#FRAGMENTS.append( [NL(self.chains[1])-H2O+PROTON, 1000, "sn2 acylium ion"])
			# 		#FRAGMENTS.append( [NL(self.chains[1])-H2O-H2O+PROTON, 1000, "sn2 acylium ion - H2O"])
	
			
			elif adduct in NEG_ADDUCTS:  
				if adduct == "[M+Cl]-" or adduct == "[M+FA]-":
					FRAGMENTS.append( [PREC-ADDUCT[adduct]-PROTON, 1000, "NL adduct"] )
				if adduct == "[M-H]-":
					FRAGMENTS.append([PREC, 1000, "Precursor"])
					FRAGMENTS.append( [PREC - H2O, 3, "NL water"] )
					FRAGMENTS.append( [PREC - MW("CO2"), 5, "NL CO2"] )
					FRAGMENTS.append( [MW("C11H19NO9")-PROTON, 15, "N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [MW("C11H19NO9")-H2O-PROTON, 388, "deoxy-N-Ac-Neuraminic Acid"] )

					FRAGMENTS.append( [PREC-MW("C11H19NO9")+H2O, 76, "NL N-Ac-Neuraminic Acid"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")+2*H2O, 28, "pre - N-Ac-Neuraminic Acid - N Ac Glucosamine"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-MW("C6H10O5")+2*H2O, 15, "pre - N-Ac-Neuraminic Acid - N Ac Glucosamine - glucose"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")+2*H2O, 33, "pre - N-Ac-Neuraminic Acid - N Ac Glucosamine - glucose - galactose"] )
					FRAGMENTS.append( [PREC-MW("C11H19NO9")-MW("C8H15NO6")-2*MW("C6H10O5")+H2O, 2, "pre - N-Ac-Neuraminic Acid - N Ac Glucosamine - glucose - galactose - water"] )

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



# calicolipidlibrary.print_spectrum("AcGM2",[[18,1,1], [18,0,0]], "[M+H]+")

# x  = AcGM2("AcGM2",[[18,1,1], [18,0,0]],  adduct="[M+H]+")
# print x.printNist()

