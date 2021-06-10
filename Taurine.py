from lipidRules import *

class Taurine(singleAcyl):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
		
				FRAGMENTS.append([PREC, 35, "precursor"])
				FRAGMENTS.append( [MW("C2H7NO3S")+ADDUCT[adduct],1000,"Taurine"])
			
		
		
			
			if adduct in NEG_ADDUCTS:
				FRAGMENTS.append([PREC, 1000, "precursor"])
		
		
				FRAGMENTS.append( [PREC-MW("C2H7NO3S"),1000,"NL Taurine"])
				FRAGMENTS.append( [MW("C2H7NO3S")-PROTON,27,"Taurine"])
				FRAGMENTS.append( [MW("C2H7NO3S")-MW("NH3")-PROTON,40,"Taurine - NH3"])
				FRAGMENTS.append( [MW("SO3H")-PROTON,45,"sulfate"])

			
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




x  = Taurine("Taurine",[[18,0,0]],  adduct="[M+H]+")
print x.printNist()

