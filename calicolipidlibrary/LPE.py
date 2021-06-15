from lipidRules import *

class LPE(LysoGPL):



		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct
			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:
            

  
				if adduct == "[M+H]+":
					FRAGMENTS.append([PREC, 23, "precursor"])
				
					FRAGMENTS.append( [PREC-H2O, 37, "NL H2O"] )
					FRAGMENTS.append( [PREC-MW("C2H5N")-H2O, 35, "NL C2NH5 + H2O"] )
					FRAGMENTS.append( [PREC-MW("C2H8NO4P"), 1000, "NL ethanolamine-P"] )
					FRAGMENTS.append( [PREC-MW("C3H9O6P"), 90, "NL glycerol-P"] )
					FRAGMENTS.append( [MW("C5H14NO6P")+ADDUCT[adduct], 25, "glycerol+phosphoethanolamine"] )
					FRAGMENTS.append( [MW("C5H14NO6P")-H2O+ADDUCT[adduct], 16, "glycerol+phosphoethanolamine - H2O"] )
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct], 10, "glycerol-P"] )
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct]-H2O, 60, "glycerol-P - H2O"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct], 13, "ethanolamine-P"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct]-H2O, 10, "ethanolamine-P - H2O"] )
					FRAGMENTS.append( [MW("H3O4P")+ADDUCT[adduct], 15, "PO4H4+"] )
					FRAGMENTS.append( [MW("C2H7NO")+ADDUCT[adduct], 234, "ethanolamine"] )
					for i in range(0,len(self.chains)):
						chain = str(i+1)
						if self.chains[i][0] == 0:
							continue

						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O-H2O, 9, "sn" + chain + " acylium ion - water"])
						FRAGMENTS.append( [NL(self.chains[i])+PROTON-H2O, 24, "sn" + chain + " acylium ion"])

				if adduct in ["[M+Na]+", "[M+K]+", "[M+Li]+", "[M+NH4]+"]:
					FRAGMENTS.append([PREC, 367, "precursor"])

					FRAGMENTS.append( [PREC-MW("NH3"), 3, "NL NH3"] )				
					FRAGMENTS.append( [PREC-H2O, 10, "NL water"] )
					FRAGMENTS.append( [MW("PO4H3")+ADDUCT[adduct], 215, "phosphoric acid + adduct"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct], 74, "ethanolamine-P + adduct"] )
					FRAGMENTS.append( [PREC-MW("C2H5N"), 1000, "NL deoxyethanolamine"] )
					FRAGMENTS.append( [PREC-MW("C2H7NO"), 240, "NL ethanolamine"] )
					
					FRAGMENTS.append( [PREC-MW("C2H8NO4P"), 12, "NL ethanolamine-P"] )
					FRAGMENTS.append( [PREC-MW("C2H8NO4P")-ADDUCT[adduct]+PROTON, 397, "NL ethanolamine-P + ADDUCT"] )

					FRAGMENTS.append( [PREC-MW("C2H6NO3P"), 44, "NL deoxyethanolamine-P"] )       	
					
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct]-H2O, 47, "ethanolamine-P - H2O"] )
					FRAGMENTS.append( [MW("C3H9O6P")+ADDUCT[adduct], 121, "glycerol-P"] )
					for i in range(0,(len(self.chains))):
						chain = str(i+1)
						if self.chains[i][0] == 0:
							continue
						FRAGMENTS.append( [NL(self.chains[i])-H2O+PROTON, 5, "sn1 acylium ion"])
						FRAGMENTS.append( [PREC-MW("C2NH5")-NL(self.chains[i]), 40, "NL sn1 + deoxyethanolamine"] ) 

			
			elif adduct in NEG_ADDUCTS:

				if self.chains[0][0] != 0:
					i = 0

					FRAGMENTS.append([PREC, 152, "precursor"])

					FRAGMENTS.append( [MW("C3H9O6P")-H2O-PROTON, 10,  "glycerol3P - water"] )
					FRAGMENTS.append( [MW("HPO3")-PROTON,  42,  "phosphite"] )
					FRAGMENTS.append( [MW("C2H8NO4P")-H2O+ADDUCT[adduct],  2,  "deoxyethanolamine-P + adduct"] )
					FRAGMENTS.append( [MW("C2H8NO4P")+ADDUCT[adduct],  32,  "ethanolamine-P + adduct"] )

					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 1000, "sn1"])
					FRAGMENTS.append( [PREC - NL(self.chains[i])+H2O, 24, "NL sn1 ketene"])
					FRAGMENTS.append( [PREC - NL(self.chains[i]), 100, "NL sn1"])

				elif self.chains[1][0] != 0:
					i = 1
					FRAGMENTS.append([PREC, 78, "precursor"])

					FRAGMENTS.append([MW("C3H9O6P") - H2O - PROTON, 2, "glycerol3P - water"])
					FRAGMENTS.append([MW("HPO3") - PROTON, 23, "phosphite"])
					#FRAGMENTS.append([MW("C2H8NO4P") - H2O + ADDUCT[adduct], 2, "deoxyethanolamine-P + adduct"])
					FRAGMENTS.append([MW("C2H8NO4P") + ADDUCT[adduct], 15, "ethanolamine-P + adduct"])

					FRAGMENTS.append([NL(self.chains[i]) - PROTON, 1000, "sn1"])
					FRAGMENTS.append([PREC - NL(self.chains[i]) + H2O, 45, "NL sn1 ketene"])
					FRAGMENTS.append([PREC - NL(self.chains[i]), 20, "NL sn1"])



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

# x  = LPE("LPE",[[18,0,0],[0,0,0]],  adduct="[M+H]+")
# print x.printNist()
#
# calicolipidlibrary.print_spectrum("LPE", [[18,0,0],[0,0,0]], "[M+H]+")