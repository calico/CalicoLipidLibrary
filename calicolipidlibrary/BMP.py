from lipidRules import *

class BMP(GPL):

		pos_adduct_set = ["[M+H]+"]

		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]


			if adduct in POS_ADDUCTS:  
				if adduct == "[M+Na]+" or adduct == "[M+K]+":
					FRAGMENTS.append([PREC, 1000, "precursor"])
					FRAGMENTS.append([PREC-H2O, 64, "NL water"])				
					FRAGMENTS.append([MW("C3H9O6P")+ADDUCT[adduct], 60, "glycerol-P + adduct"])							
					FRAGMENTS.append([MW("C3H9O6P")-H2O+ADDUCT[adduct], 54, "deoxyglycerol-P + Adduct"])							
					FRAGMENTS.append([MW("PO4H3")+ADDUCT[adduct], 36, "phosphoric acid + Adduct"])							


					for i in range(0,len(self.chains)):
						chain = str(i + 1)						
						FRAGMENTS.append([PREC-NL(self.chains[i])-MW("C3H4O"), 378, "NL (sn"+chain+" + deoxyglycerol)"])
						FRAGMENTS.append([PREC-NL(self.chains[i])-MW("C3H6O2"), 110, "NL (sn"+chain+" + glycerol)"])		
						FRAGMENTS.append([PREC-NL(self.chains[i])-MW("C3H4O")-MW("PO4H3")-ADDUCT[adduct]+PROTON, 188, "NL (sn"+chain+" + deoxyglycerol + adduct + phosphate)"])		
						FRAGMENTS.append([PREC-NL(self.chains[i])-MW("C3H4O")-MW("PO4H3"), 5, "NL (sn"+chain+" + deoxyglycerol + phosphate)"])							
						FRAGMENTS.append([PREC-NL(self.chains[i])-MW("C3H4O")-MW("PO4H3") + H2O, 18, "NL (sn"+chain+" + deoxyglycerol + phosphate) + H2O"])							
	
						FRAGMENTS.append([PREC-NL(self.chains[i]), 2, "NL sn" + chain])				
						FRAGMENTS.append([PREC-NL(self.chains[i])+H2O, 17, "NL sn"+chain+" ketene"])				


				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 1, "precursor"])
					FRAGMENTS.append([PREC-H2O, 1, "NL water"])				
					FRAGMENTS.append([PREC-2*H2O, 158, "NL - 2x water"])				
					FRAGMENTS.append([PREC - MW("C3H9O6P"), 27, "NL glycerol-P"])							



					for i in range(0,len(self.chains)):
						chain = str(i + 1)	
						
						FRAGMENTS.append([NL(self.chains[i])+MW("C3H4O")+PROTON, 500, "(sn"+chain+" + deoxyglycerol)"])
						FRAGMENTS.append([PREC-NL(self.chains[i])-MW("C3H6O2"), 2, "NL (sn"+chain+" + glycerol)"])		
						FRAGMENTS.append([PREC-NL(self.chains[i])-H2O, 3, "NL (sn" + chain +" ketene)"])				





			
			
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 501, "precursor"])

				FRAGMENTS.append( [MW("C3H7O5P")-PROTON, 41, "C3H6O5P-"] )

				FRAGMENTS.append( [MW("C6H13O7P")-PROTON, 5, "diglycerol-phosphate"] )
				FRAGMENTS.append( [MW("C6H11O6P")-PROTON, 1, "diglycerol-phosphate - water"] )

				for i in range(0,len(self.chains)):
					chain = str(i + 1);		
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 500, "sn"+chain])
					FRAGMENTS.append( [PREC-NL(self.chains[i]), 12, "NL (sn"+chain])
					FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O, 3, "NL sn"+chain+" ketene"])
					FRAGMENTS.append( [PREC-NL(self.chains[i])+H2O-MW("C3H8O3"), 7, "NL (sn"+chain+" ketene + glycerol)"])
					FRAGMENTS.append( [PREC-NL(self.chains[i])+2*H2O-MW("C3H8O3"), 1, "NL (sn"+chain+" ketene + deoxyglycerol)"])

								

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



x  = BMP("BMP",[[18,1,0], [18,1,0]],  adduct="[M+H]+")
print x.printNist()



x  = BMP("BMP",[[10,0,0], [12,4,0]],  adduct="[M-H]-")
print x.printNist()

