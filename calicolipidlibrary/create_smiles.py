

HEADGROUP_SMILES = {
    # backbone,   #linkage
    "PC": "(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)",
    "PE": "(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)"
}

ACYL_SMILES = {}

for c in range(2,35):
    for d in range (0,8):
        for h in range(0,1): #this is just zero for the time being.
            #add checks to continue of dbs are too high
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            smiles_str = "C" * c
            match d:
                case 1:
                    if c < 13:
                        smiles_str = smiles_str[0:-8] + "/C=C\\" + smiles_str[-6:]
                    elif:
                        smiles_str = smiles_str[0:4] + "/C=C\\" + smiles_str[6:]



            ACYL_SMILES.update({["acyl", c, d, h]: smiles_str})

# [[type, C, DB, H], "SMILE STRING]




}