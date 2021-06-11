from lipidRules import *

class Carn(singleAcyl):


	

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]


			if adduct in POS_ADDUCTS:

				FRAGMENTS.append([PREC, 837, "precursor"])
				FRAGMENTS.append( [MW("C3H9N")+PROTON,83,"trimethylamine"])
			
				FRAGMENTS.append( [PREC-MW("C3H9N"),129,"NL trimethylamine"])
		
				FRAGMENTS.append( [MW("C4H4O2")+ADDUCT[adduct],890,"C4H5O2 + adduct "])
				FRAGMENTS.append( [MW("C7H13NO2")+ADDUCT[adduct],32,"deoxycarnitine"])                 
				FRAGMENTS.append( [NL(self.chains[0])-H2O+ADDUCT[adduct],40,"acylium ion "])                 
			


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


x  = Carn("Carn",[[12,0,0]],  adduct="[M+H]+")
print x.printNist()



x  = Carn("Carn",[[16,0,0]],  adduct="[M+H]+")
print x.printNist()

