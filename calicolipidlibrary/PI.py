from lipidRules import *

class PI(GPL):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct
		
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
				FRAGMENTS.append([PREC, 553, "precursor"])
	   
				FRAGMENTS.append( [PREC-MW("C6H13O9P")-ADDUCT[adduct]+PROTON, 134, "NL inositol-P + adduct "] )

				if adduct == "[M+H]+":
					FRAGMENTS.append( [NL(self.chains[0])-H2O+PROTON, 7, "sn1 acylium ion"])
					FRAGMENTS.append( [NL(self.chains[1])-H2O+PROTON, 7, "sn2 acylium ion"])



				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append( [PREC-MW("C6H13O9P"), 161, "NL inositol-phosphate"] )
					FRAGMENTS.append( [PREC-MW("C6H13O9P")+H2O, 15, "NL deoxyinositol "] )      
					FRAGMENTS.append( [MW("C6H13O9P")+ADDUCT[adduct], 1000, "inositol-phosphate + adduct"] )        
					FRAGMENTS.append( [MW("C6H13O9P")-H2O+ADDUCT[adduct], 37, "inositol-phosphate - water + adduct"] )        
	
					i = 0
					chain = str(i+1)
					FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 5, "sn" + chain + " acylium ion"])
					FRAGMENTS.append( [PREC - MW("C6H13O9P") - NL(self.chains[i])+H2O-ADDUCT[adduct]+PROTON, 3, "NL sn"+chain+" + inositol-P + adduct"])
					FRAGMENTS.append( [PREC - MW("C6H12O6") - NL(self.chains[i])+H2O, 7, "NL (sn"+chain+" + inositol-P)"])

					i = 1
					chain = str(i+1)
					FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 10, "sn" + chain + " acylium ion"])
					FRAGMENTS.append( [PREC - MW("C6H13O9P") - NL(self.chains[i])+H2O-ADDUCT[adduct]+PROTON, 37, "NL sn"+chain+" + inositol-P + Adduct"])
					FRAGMENTS.append( [PREC - MW("C6H12O6") - NL(self.chains[i])+H2O, 14, "NL (sn"+chain+" + inositol-P)"])


		
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])

		
				FRAGMENTS.append( [MW("C3H9O6P")-H2O+ADDUCT[adduct], 66,  "glycerol3P-H2O"] ) 
				FRAGMENTS.append( [MW("C3H9O6P")-H2O+ADDUCT[adduct], 2,  "glycerol3P - 2xH2O"] ) 				
				FRAGMENTS.append( [MW("H3PO4")-PROTON,  8,  "H2PO4-"] )  #sometimes this mass is too low

				FRAGMENTS.append( [MW("C9H19O11P")-H2O-PROTON,  213,  "glycerol-inositol-phosphate - H2O"] )
				FRAGMENTS.append( [MW("C9H19O11P")-H2O-H2O-PROTON,  46,  "glycerol-inositol-phosphate - 2xH2O"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON,  21,  "inositol-phosphate"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O, 213,  "inositol-phosphate - H2O"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O-H2O,  46,  "inositol-phosphate - 2xH2O"] )

				i = 0
				chain = str(i+1)
				if (self.chains[i][1] > 3): #add decarboxylation if PUFA
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 280, "sn"+chain+ " fatty acyl"])
					FRAGMENTS.append( [NL(self.chains[i])-PROTON-MW("CO2"), 30, "sn"+chain+ " fatty acyl - CO2"])
				else: 
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 309, "sn"+chain+ " fatty acyl"])

				FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 4, "NL sn"+chain])
				FRAGMENTS.append( [PREC - NL(self.chains[i]), 23, "NL sn"+chain+ "-H2O"])
				FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O - MW("C6H12O6") +H2O ,1, "NL sn"+chain+ "- H2O - inositol"])
				FRAGMENTS.append( [PREC - NL(self.chains[i]) - MW("C6H12O6") +H2O, 19, "NL sn"+chain+ "- 2xH2O - inositol"])

				i = 1
				chain = str(i+1)
				if (self.chains[i][1] > 3): #add decarboxylation if PUFA
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 264, "sn"+chain+ " fatty acyl"])
					FRAGMENTS.append( [NL(self.chains[i])-PROTON-MW("CO2"), 25, "sn"+chain+ " fatty acyl - CO2"])
				else: 
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 293, "sn"+chain+ " fatty acyl"])
				FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 16, "NL sn"+chain])
				FRAGMENTS.append( [PREC - NL(self.chains[i]), 63, "NL sn"+chain+ "-H2O"])
				FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O - MW("C6H12O6") +H2O ,5, "NL sn"+chain+ "-H2O-INOS"])
				FRAGMENTS.append( [PREC - NL(self.chains[i]) - MW("C6H12O6") +H2O, 65, "NL sn"+chain+ "-H2O-H2O-INOS"])

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




# x  = PI("PI",[[20,4,0], [16,0,0]],  adduct="[M+Na]+")
# print x.printNist()
# x  = PI("PI",[[16,0,0], [20,4,0]],  adduct="[M-H]-")
# print x.printNist()

# calicolipidlibrary.print_spectrum("PI", [[20,4,0], [16,0,0]], "[M+Na]+")