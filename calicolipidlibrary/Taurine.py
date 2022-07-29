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



# x  = Taurine("Taurine",[[18,0,0]],  adduct="[M+H]+")
# print x.printNist()

# calicolipidlibrary.print_spectrum("Taurine", [[18,0,0]], "[M+H]+")