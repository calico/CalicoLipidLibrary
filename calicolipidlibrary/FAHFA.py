from lipidRules import *


#assumptions:  sn2 has the free carboxyl, sn1 is the esterified FA (this matches LipidMaps nomenclature).
#For our purposes a standard FAHFA has no hydroxyls, as the esterified one is counted as part of SN1

class FAHFA(GPL):
		pos_adduct_set = ["[M+Na]+", "[M+K]+", "[M+NH4]+"]

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in NEG_ADDUCTS:
				FRAGMENTS.append([PREC, 81, "precursor"])

				FRAGMENTS.append( [NL(self.chains[1])-MW("H2") - PROTON,256,"sn2"])
				FRAGMENTS.append( [NL(self.chains[1])-MW("H2")+H2O - PROTON,188,"sn2 + water"])
				FRAGMENTS.append( [NL(self.chains[0]) - PROTON,1000,"sn1"])


			if adduct in POS_ADDUCTS:
				if adduct == "[M+Na]+" or adduct == "[M+K]+" or adduct == "[M+NH4]+":

					FRAGMENTS.append([PREC, 1000, "precursor"])

					FRAGMENTS.append( [NL(self.chains[1])-MW("H2") +ADDUCT[adduct],21,"sn2 + adduct"])
					FRAGMENTS.append( [NL(self.chains[1])-MW("H2")+H2O +ADDUCT[adduct],1,"sn2 + water + adduct"])
					FRAGMENTS.append( [NL(self.chains[0]) +ADDUCT[adduct],8,"sn1 + adduct"])
		
			


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




# x  = FAHFA("FAHFA",[[18,0,0], [18,1,0]],  adduct="[M+Na]+")
# print x.printNist()
#
# x  = FAHFA("FAHFA",[[16,0,0], [18,0,0]],  adduct="[M+Na]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("FAHFA", [[18,0,0], [18,1,0]], "[M+Na]+")