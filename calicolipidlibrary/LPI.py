from lipidRules import *

class LPI(LysoGPL):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				# no sign of M+H from standards
				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 2, "precursor"])


					FRAGMENTS.append( [PREC-H2O, 17, "NL water"])
					FRAGMENTS.append( [PREC-2*H2O, 44, "NL 2x water"])

					FRAGMENTS.append( [PREC-MW("C6H12O6"), 26, "NL inositol"] )
					FRAGMENTS.append( [PREC-MW("C6H13O9P"), 1000, "NL inositol phosphate"] )
					FRAGMENTS.append( [MW("C6H13O9P")+ADDUCT[adduct], 22, "inositol-phosphate + adduct"] )
					FRAGMENTS.append( [MW("C6H13O9P")-H2O+ADDUCT[adduct], 4, "inositol-phosphate - water"] )
					FRAGMENTS.append( [MW("C6H12O6")-H2O+ADDUCT[adduct], 3, "inositol - water"] )
					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 36, "phosphoric acid"] )
					FRAGMENTS.append( [MW("C3H9O6P")-H2O+PROTON, 132, "deoxyglycerol-P"] )
					FRAGMENTS.append( [MW("C3H9O6P")+PROTON, 16, "glycerol-P"] )
			
					FRAGMENTS.append( [NL(self.chains[0])-H2O+PROTON, 38, "sn1 acylium ion"])
					FRAGMENTS.append( [PREC - NL(self.chains[0])+H2O, 27, "NL ketene"])
				 

	
				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 1000, "precursor"])
				
					FRAGMENTS.append( [PREC-H2O, 38, "NL water"])
					FRAGMENTS.append( [PREC-MW("C6H12O6"), 335, "NL inositol"] )      
					FRAGMENTS.append( [PREC-MW("C6H10O5"), 75, "NL deoxy-inositol"] )  				
					FRAGMENTS.append( [PREC-MW("C6H13O9P"), 3, "NL inositol phsophate"] )
					FRAGMENTS.append( [PREC-MW("C6H13O9P")+H2O, 30, "NL inositol-phosphate + water "] )      
					FRAGMENTS.append( [MW("C6H13O9P")+ADDUCT[adduct], 377, "inositol-phosphate + adduct"] )
					FRAGMENTS.append( [MW("C6H13O9P")-H2O+ADDUCT[adduct], 134, "inositol-phosphate - water + adduct"] )
					FRAGMENTS.append( [MW("C6H12O6")+ADDUCT[adduct], 55, "inositol + adduct"] )					
					FRAGMENTS.append( [MW("C6H12O6")-H2O+ADDUCT[adduct], 5, "inositol - water + adduct"] )					
					FRAGMENTS.append( [MW("H3PO4")+ADDUCT[adduct], 5, "phosphoric acid + adduct"] )
					FRAGMENTS.append( [MW("C3H9O6P")-H2O+ADDUCT[adduct], 53, "glycerol-P - water + adduct"] )
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct], 14, "glycerol-P  + adduct"] )
					
					# additional peak at 176.9920171 seen in two standards, but cannot rationalize, ion count = 53
					for i in range(0,len(self.chains)):
						chain = str(i+1)
						if self.chains[i][0] == 0:
							continue
						FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 5, "sn1 acylium ion"])
						FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 27, "NL sn1 ketene"])
				 
			
			elif adduct in NEG_ADDUCTS:

				if self.chains[0][0] != 0:
					i = 0
				else:
					i = 1


				FRAGMENTS.append([PREC, 1000, "precursor"])

				FRAGMENTS.append( [MW("C3H9O6P")-H2O+ADDUCT[adduct], 291,  "glycerol3P - water"] )
				FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct], 5,  "glycerol3P"] )
				FRAGMENTS.append( [MW("H3PO4")-PROTON,  10,  "phosphate"] )
				FRAGMENTS.append( [MW("HPO3")-PROTON,  34,  "phosphite"] )
				FRAGMENTS.append( [MW("C6H12O6")-PROTON,  1,  "inositol"] )


				FRAGMENTS.append( [MW("C9H19O11P")-H2O-PROTON,  137,  "glycerol-phosphoinositol - H2O"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON,  7,  "inositol-phosphate"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O,  228,  "inositol-phosphate - H2O"] )
				FRAGMENTS.append( [MW("C6H13O9P")-PROTON-H2O-H2O,  13,  "inositol-phosphate - 2xH2O"] )
				FRAGMENTS.append( [PREC-MW("C6H10O5"),  5,  "NL deoxyinositol"] )
				FRAGMENTS.append( [PREC-MW("C6H12O6"),  51,  "NL inositol"] )


				chain = str(i+1)
				FRAGMENTS.append( [NL(self.chains[i])-PROTON, 505, "sn"+chain])
				FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 5, "NL sn"+chain])


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




x  = LPI("LPI",[[16,0,0],[0,0,0]],  adduct="[M+Na]+")
print x.printNist()
x  = LPI("LPI",[[20,4,0],[0,0,0]],  adduct="[M-H]-")
print x.printNist()

