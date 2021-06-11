from lipidRules import *

class CDP_DG(GPL):

		pos_adduct_set = ["[M+H]+", "[M+Na]+", "[M+K]+"]


		def theoreticalDigest(self):
			FRAGMENTS = []
			adduct = self.adduct

			
			MASS = MW_list(self.MF())
			PREC  = MASS + ADDUCT[self.adduct]

			if adduct in POS_ADDUCTS:

				if adduct == "[M+H]+":      
					FRAGMENTS.append([PREC, 24, "precursor"])
					FRAGMENTS.append([MW("C4H5N3O")+PROTON,1000, "cytosine"])
					FRAGMENTS.append([MW("C9H13N3O5")-2*H2O+PROTON, 81,"deoxyCytidine - water "])


					FRAGMENTS.append([PREC - MW("C9H15N3O11P2"), 703,"NL deoxyCDP"])
					FRAGMENTS.append([MW("C9H15N3O11P2")+PROTON, 47,"deoxyCDP"])
					for i in range(0,len(self.chains)):
						chain = str(i+1)
						FRAGMENTS.append([NL(self.chains[i])-H2O+PROTON, 14, "sn"+chain+" acylium"])
					
						

				if adduct == "[M+Na]+" or adduct == "[M+K]+":  
					FRAGMENTS.append([PREC, 149, "precursor"])
					FRAGMENTS.append([MW("C4H5N3O")+ADDUCT[adduct],140, "cytosine + adduct"])
					FRAGMENTS.append([MW("C9H13N3O5")-H2O+PROTON, 81,"deoxycytidine"])
					FRAGMENTS.append([MW("P2O7H4")+ADDUCT[adduct], 325,"pyrophosphate + Adduct"])
					FRAGMENTS.append([MW("C5H12O11P2")+ADDUCT[adduct]-2*H2O, 289,"ribose-pyrophosphate - 2x water"])
					FRAGMENTS.append([MW("C5H12O11P2")+ADDUCT[adduct]-H2O, 9,"ribose-pyrophosphate -  water"])
					FRAGMENTS.append([MW("C5H11O8P")+ADDUCT[adduct]-2*H2O, 114,"ribose-phosphate -  2x water"])
					FRAGMENTS.append([MW("C5H11O8P")+ADDUCT[adduct]-3*H2O, 30,"ribose-phosphate -  3x water"])



					FRAGMENTS.append([PREC - MW("C9H15N3O11P2")-ADDUCT[adduct]+PROTON, 304,"NL deoxyCDP + adduct"])
					FRAGMENTS.append([PREC - MW("C9H13N3O5")+H2O, 25,"NL deoxyCytidine"])
					
					FRAGMENTS.append([MW("C9H15N3O11P2")+ADDUCT[adduct], 1000,"deoxyCDP + adduct"])
					FRAGMENTS.append([MW("C9H14N3O8P")+ADDUCT[adduct], 157,"deoxyCMP + adduct"])
					FRAGMENTS.append([MW("C9H14N3O8P")-H2O+ADDUCT[adduct], 99,"deoxyCMP - water + adduct"])
					
					
					for i in range(0,len(self.chains)):
						chain = str(i+1)
						FRAGMENTS.append([NL(self.chains[i])-H2O+PROTON, 11, "sn"+chain+" acylium"])
						FRAGMENTS.append([NL(self.chains[i])+MW("C3H4O")+PROTON, 4, "sn"+chain+"+ deoxyglycerol"])
						
						
						
						
			
			elif adduct in NEG_ADDUCTS:  
				FRAGMENTS.append([PREC, 1000, "precursor"])
				FRAGMENTS.append([MW("C9H15N3O11P2")-H2O-PROTON, 238,"CDP - water"])
				FRAGMENTS.append([MW("C9H14N3O8P")-H2O-PROTON, 69,"CMP - water"])
				FRAGMENTS.append([MW("C9H14N3O8P")-PROTON, 12,"CMP "])


				#big peak @ 272.9522743, not depentend on acyl groups
				
				FRAGMENTS.append([PREC - MW("C9H13N3O5"), 482,"NL Cytidine"])
				FRAGMENTS.append([MW("P2O7H4")-PROTON, 7,"pyrophosphate"])
				FRAGMENTS.append([MW("P2O7H4")-H2O-PROTON, 113,"pyrophosphate - water"])


				for i in range(0,len(self.chains)):
					chain = str(i+1)
					FRAGMENTS.append( [NL(self.chains[i])-PROTON, 16, "sn"+chain] )								
					FRAGMENTS.append( [PREC - NL(self.chains[i]), 8, "NL sn"+chain] )								
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - MW("C9H13N3O5"), 18, "NL sn"+chain +" + Cytidine"])								
					FRAGMENTS.append( [PREC - NL(self.chains[i]) - MW("C9H13N3O5")+H2O, 6, "NL sn"+chain +" + deoxycytidine"])								


			   

				

			
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




x  = CDP_DG("CDP_DG",[[18,1,0], [18,1,0]],  adduct="[M+Na]+")
print x.printNist()

x  = CDP_DG("CDP_DG",[[16,0,0], [16,0,0]],  adduct="[M+Na]+")
print x.printNist()

