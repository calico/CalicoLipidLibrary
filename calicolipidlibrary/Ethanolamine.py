from lipidRules import *

class Ethanolamine(singleAcyl):


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]
			FRAGMENTS.append([PREC, 1000, "precursor"])

			if adduct in POS_ADDUCTS: #intensites from literature, not spectra
			
				FRAGMENTS.append( [PREC-H2O,300,"NL water"])

				FRAGMENTS.append( [PREC-MW("C2H7NO"),300,"NL ethanolamine"])
				FRAGMENTS.append( [MW("C2H7NO")+PROTON,1000,"Ethanolamine"])
				
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




# x  = Ethanolamine("Ethanolamine",[[16,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("Ethanolamine", [[16,0,0]], "[M+H]+")