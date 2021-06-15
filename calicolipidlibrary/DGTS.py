from lipidRules import *

class DGTS(GPL):
		pos_adduct_set = ["[M+H]+"]

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
			
			
				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append([MW("C7H15O3N")+PROTON, 11, "trimethylhomoserine"])				
					FRAGMENTS.append([MW("C7H15O3N")-H2O+PROTON, 35, "deoxytrimethylhomoserine"])				
					FRAGMENTS.append([MW("C7H15O3N")+MW("C3H8O3")-H2O+PROTON, 247, "glycerotrimethylhomoserine"])				
					FRAGMENTS.append([MW("C7H15O3N")+MW("C3H8O3")-2*H2O+PROTON, 3, "glycerotrimethylhomoserine - H2O"])				
					FRAGMENTS.append([MW("C7H15O3N")+MW("C3H8O3")-3*H2O+PROTON, 6, "glycerotrimethylhomoserine - 2*H2O"])				


					for i in range(0,len(self.chains)):
						chain = str(i+1)

						FRAGMENTS.append( [PREC - NL(self.chains[i])-H2O, 210, "NL sn" + chain +" acylium ion"])
						FRAGMENTS.append( [PREC - NL(self.chains[i]), 66, "NL sn" + chain +" acylium ion"])



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




# x  = DGTS("DGTS",[[16,0,0], [16,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("DGTS", [[16,0,0], [16,0,0]], "[M+H]+")