from lipidRules import *

#NEEDS WORK....how to generate acyl combinations such that we get the right fragments....
#for any 3 FA combination we need to cover the basis where each individual FA is on the glycerol by itself

class LysoCL(LysoCardioLipin):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]
			FRAGMENTS.append([PREC, 1000, "precursor"])

			
			
			if adduct in NEG_ADDUCTS:  
				FRAGMENTS.append ( [PREC-ADDUCT[adduct]-PROTON-H2O,1000,"NL water"])
				if adduct == "[M-H]-":
					FRAGMENTS.append( [MW("C3H9O6P")-PROTON, 1000, "Glycerol-P"] )
					FRAGMENTS.append( [MW("C3H9O6P")-H2O-PROTON, 1000, "Glycerol-P - water"] )					
					for i in range(0,len(self.chains)):
						chain = str(i + 1);

						FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn"+chain])
						FRAGMENTS.append( [PREC-NL(self.chains[i]), 1000, "nL (sn"+chain+" + water)"])
						FRAGMENTS.append( [NL(self.chains[i])+MW("C3H7O5P")-PROTON, 1000, "sn"+chain+" LPA"])
						FRAGMENTS.append( [NL(self.chains[i])-H2O+MW("C3H7O5P")-PROTON, 1000, "sn"+chain+" LPA - water"])
						

					FRAGMENTS.append( [PREC-NL(self.chains[0])-NL(self.chains[1])-MW("C3H9O6P"), 1000, "NL glycerol-P + sn1+sn2" ])
					FRAGMENTS.append( [PREC-NL(self.chains[2])-MW("C3H9O6P"), 1000, "NL glycerol-P + sn3+sn4" ])

					FRAGMENTS.append( [NL(self.chains[0])+NL(self.chains[1])+MW("C3H9O6P")-PROTON, 1000, "glycerol-P + sn1+sn2" ])
					FRAGMENTS.append( [NL(self.chains[2])++MW("C3H9O6P")-PROTON, 1000, "glycerol-P + sn3+sn4" ])
		
					FRAGMENTS.append( [PREC-NL(self.chains[0])-NL(self.chains[1])-MW("C3H8O3"), 1000, "NL glycerol + sn1+sn2" ])
					FRAGMENTS.append( [PREC-NL(self.chains[2])-MW("C3H8O3"), 1000, "NL glycerol + sn3+sn4" ])

				
			return(FRAGMENTS)

		def generateLibrary(self, target=None, mode="pos"):
			if target: handle = open(target, 'a+')
			if mode == "pos":
				adduct_set = self.pos_adduct_set
			elif mode == "neg":
				adduct_set = self.neg_adduct_set
			# parent = self.__bases__[0]
			class_name = self.__class__.__name__
			for c in self.chain_sets:
				for adduct in adduct_set:
					self.set_chains_and_adduct(class_name, c[:3], adduct=adduct)
					content = self.printNist()
					if target:
						handle.write(content)
					else:
						sys.stdout.write(content)
			if target: handle.close()




# x  = LysoCL("LysoCL",[[16,0,0], [18,1,0], [18,2,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("LysoCL", [[16,0,0], [18,1,0], [18,2,0]], "[M+H]+")