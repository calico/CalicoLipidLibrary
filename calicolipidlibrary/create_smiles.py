from rdkit import Chem
import copy
import numpy as np

#this create default positions for dbs.  specifics modified after
#first dimenstion in number of carbons, second is number of double bonds
acyl_db_table = {

}

for c in range(10,37):
    acyl_db_table[c] = {}
    for d in range(0,7):
        if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
            continue
        match d:  # #db
            case 0:
                acyl_db_table[c][d] = []
            case 1:
                if c < 10:
                    acyl_db_table[c][d] = [4]
                else:
                    acyl_db_table[c][d] = [9]
            case 2:
                if c < 14:
                    acyl_db_table[c][d] = [4,7]
                else:
                    acyl_db_table[c][d] = [9 ,12]
            case 3:
                acyl_db_table[c][d] = [5, 9 ,12]
            case 4:
                acyl_db_table[c][d] = [5,8,11,14]
            case 5:
                acyl_db_table[c][d] = [5,8,11,14,17]
            case 6:
                acyl_db_table[c][d] = [4,7,10,13,16,19]



#make sure db placement is in line with most abundant isomer of a specific acyl (overwrites above)

acyl_db_table[20][1] = [11]
acyl_db_table[22][1] = [13]
acyl_db_table[24][1] = [15]
acyl_db_table[20][3] = [8,11,14]
acyl_db_table[18][4] = [6,9,12,15]
acyl_db_table[22][4] = [7,10,13,16]
acyl_db_table[24][5] = [9,12,15,18,21]
acyl_db_table[24][6] = [6,9,12,15,18,21]

alkyl_db_table = copy.deepcopy(acyl_db_table)
#modify so first db is always at 1, to make plasmalogens
#for 2+ db, make plasmalogen, plus location for acyl double bonds for on less db
for c in range(10,37):
    for d in range(1,7):
        if d == 1:
            alkyl_db_table[c][d] = [1]
        else:
            alkyl_db_table[c][d] =[1] + alkyl_db_table[c][d-1]


sbase_db_table = copy.deepcopy(acyl_db_table)

#deviations in double bond placement in LCB bases on sphingolipids
sbase_db_table[16][1] = [1]
sbase_db_table[17][1] = [1]
sbase_db_table[18][1] = [1]
sbase_db_table[19][1] = [1]
sbase_db_table[20][1] = [1]
sbase_db_table[18][2] = [1,10]
#add others as seen in real data


#take in any smiles and return canonical smiles
def convert_to_canonical_SMILES(smiles_str):
    can_smiles_str = Chem.rdmolfiles.MolToSmiles(Chem.MolFromSmiles(smiles_str))
    return can_smiles_str

def add_dbs(smile_str, position_list, type): #positions use delta-# notation, so distance from the carboxylic acid
    smile_str_list = list(smile_str)
    if len(position_list) == 0:
        return smile_str
    for db in position_list:
        if type == "trans":
            smile_str_list[-(db+1)] = "/C=C/"
        elif type == "cis":
            smile_str_list[-(db+1)] = "/C=C\\"
        smile_str_list[-db] = ""
    new_smile_str = ''.join(smile_str_list)
    return new_smile_str


def get_TG_SMILES(lipidclass, sn_SMILES):
    sn1_SMILES = sn_SMILES[0]
    sn2_SMILES = sn_SMILES[1][::-1]  # reversed because at the other end of the SMILES
    sn3_SMILES = sn_SMILES[2][::-1]  # reversed because at the other end of the SMILES
    COMPLETE_SMILES = "O=C(OCC(OC(=O)" + sn1_SMILES + ")COC(=O)" + sn2_SMILES + ")" +  sn3_SMILES
    return COMPLETE_SMILES

def get_CL_SMILES(lipidclass, sn_SMILES):
    sn1_SMILES = sn_SMILES[0]
    sn2_SMILES = sn_SMILES[1][::-1]  # reversed because at the other end of the SMILES
    sn3_SMILES = sn_SMILES[2][::-1]  # reversed because at the other end of the SMILES
    sn4_SMILES = sn_SMILES[3][::-1]  # reversed because at the other end of the SMILES
    COMPLETE_SMILES = sn1_SMILES + "(=O)OC[C@@H](OC(=O)" + sn2_SMILES + ")COP(=O)(OCC(O)COP(=O)(OC[C@H](OC(=O)" + \
                        sn3_SMILES + ")COC(=O)" + sn4_SMILES + ")O)O"
    return COMPLETE_SMILES


def get_LysoCL_SMILES(lipidclass, sn_SMILES):
    sn1_SMILES = sn_SMILES[0]
    sn2_SMILES = sn_SMILES[1][::-1]  # reversed because at the other end of the SMILES
    sn3_SMILES = sn_SMILES[2][::-1]  # reversed because at the other end of the SMILES

    COMPLETE_SMILES = sn1_SMILES + "(=O)OC[C@H](COP(=O)([O-])OCC(COP(=O)([O-])OC[C@@H](COC(=O)" + \
                        sn2_SMILES + ")OC(=O)" + sn3_SMILES + ")O)O"
    return COMPLETE_SMILES


def get_BDP_SMILES(lipidclass, sn_SMILES):
    sn1_SMILES = sn_SMILES[0]
    sn2_SMILES = sn_SMILES[1][::-1]  # reversed because at the other end of the SMILES
    sn3_SMILES = sn_SMILES[2][::-1]  # reversed because at the other end of the SMILES
    sn4_SMILES = sn_SMILES[3][::-1]  # reversed because at the other end of the SMILES
    COMPLETE_SMILES = sn1_SMILES + "(=O)OC[C@@H](COP(=O)(O)OC[C@H](COC(=O)" + sn2_SMILES + ")OC(=O)" + \
                      sn3_SMILES + ")OC(=O)" + sn4_SMILES
    return COMPLETE_SMILES


def get_HemiBMP_SMILES(lipidclass, sn_SMILES):
    sn1_SMILES = sn_SMILES[0]
    sn2_SMILES = sn_SMILES[1][::-1]  # reversed because at the other end of the SMILES
    sn3_SMILES = sn_SMILES[2][::-1]  # reversed because at the other end of the SMILES
    COMPLETE_SMILES = sn1_SMILES + "(=O)OC[C@H](O)COP(=O)(O)OC[C@@H](COC(=O)" + sn2_SMILES + ")OC(=O)" + \
                      sn3_SMILES
    return COMPLETE_SMILES



def get_sphingolipid_SMILES(lipidclass, sn_SMILES):
    headgroup_SMILES = {
        # sphingolipids
        "AcGM2": "[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3NC(C)=O)[C@H](O[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)[C@H]([C@H](O)[C@H](O)CO)O3)[C@H]2O)[C@H](O)[C@H]1O)NC(=O)",
        "AcGM3": "[C@@H](O)[C@H](CO[C@@H]1OC(CO)[C@@H](O[C@@H]2OC(CO)[C@H](O)[C@H](O[C@]3(C(=O)O)CC(O)[C@@H](NC(C)=O)C([C@H](O)[C@H](O)CO)O3)C2O)[C@H](O)C1O)NC(=O)",
        "Ceramide": "[C@@H](O)[C@H](CO)NC(=O)",
        "Ceramide_P": "[C@@H](O)[C@H](COP(=O)(O)O)NC(=O)",
        # GcGM2
        "GcGM3": "[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(=O)CO)[C@H]([C@H](O)[C@H](O)CO)O3)[C@H]2O)[C@H](O)[C@H]1O)NC(=O)",
        "HexCer": "C(O)C(CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)",
        "SM": "[C@@H](O)[C@H](COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)",
        "LCB": "[C@@H](O)[C@@H](N)CO",
        "LCB_P": "[C@H]([C@H](COP(=O)(O)O)N)O",
        "Sulfatide": "[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)C(OS(=O)(=O)O)C1O)NC(=O)",
}

    if lipidclass not in headgroup_SMILES:
        return "unknown lipidclass"

    sn1_SMILES = sn_SMILES[0]
    sn2_SMILES = sn_SMILES[1][::-1]  # reversed because at the other end of the SMILES
    headgroup = headgroup_SMILES[lipidclass]

    COMPLETE_SMILES = sn1_SMILES[3:] + headgroup + sn2_SMILES[1:]

    return COMPLETE_SMILES

#currently does not work for Lysolipids with acyl at SN2 position
def get_glycerolipid_SMILES(lipidclass, sn_SMILES):

    # Headgroups here are different from that for fragmentation.
    # They generally go from the first non-carbon atom to the last non-carbon atom, or parenthesis/bracket
    headgroup_SMILES = {
        # lipidclass:   #headgroup smiles
        #phospholipids
        "Alkyl_LPC": "[C@H](COP(=O)([O-])OCC[N+](C)(C)C)O",
        "Alkyl_LPE": "[C@H](COP(=O)(O)OCCN)O",
        "Alkyl_LPS": "[C@H](COP(=O)(O)OC[C@@H](C(=O)O)N)O",
        "Alkyl_PE": "OC[C@H](COP(=O)(O)OCCN)OC(=O)",
        "Alkyl_PC": "OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)",
        "BMP": "(=O)OC[C@@H](O)COP(=O)(O)OC[C@@H](O)COC(=O)",
        "DG": "(=O)OC[C@H](CO)OC(=O)",
        "DMPE": "(=O)OC[C@H](COP(=O)(O)OCCN(C)C)OC(=O)",
        "LPA": "(=O)OCC(COP(=O)(O)O)O",
        "LPC": "(=O)OC[C@@H](COP(=O)([O-])OCC[N+](C)(C)C)O",
        "LPE": "(=O)OC[C@H](COP(=O)(O)OCCN)O",
        "LPG": "(=O)OC[C@@H](COP(=O)([O-])OC[C@H](CO)O)O",
        "LPS": "(=O)OC[C@H](COP(=O)(O)OC[C@@H](C(=O)O)N)O",
        "MG": "(=O)OCC(CO)O",
        "MMPE": "(=O)OC[C@H](COP(=O)(O)OCCNC)OC(=O)",
        "PA": "(=O)OCC(COP(=O)(O)O)OC(=O)",
        "PC": "(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)",
        "PE": "(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)",
        "PG": "(=O)OC[C@H](COP(=O)(O)OC[C@H](CO)O)OC(=O)",
        "PS": "(=O)OC[C@H](COP(=O)(O)OC[C@@H](C(=O)O)N)OC(=O)",

        #the rest
        "Carn": "(=O)O[C@H](CC(=O)[O-])C[N+](C)(C)C",
        "CE": "(=O)O[C@H]1CC[C@@]2([C@H]3CC[C@]4([C@H]([C@@H]3CC=C2C1)CC[C@@H]4[C@H](C)CCCC(C)C)C)C",
        "Ethanolamine": "(=O)NCCO",
        "FA": "(=O)O",
        "Taurine": "(=O)NCCS(=O)(=O)O",
    }

    if lipidclass not in headgroup_SMILES:
        return "unknown lipidclass"

    sn1_SMILES = sn_SMILES[0]
    sn2_SMILES = sn_SMILES[1][::-1]  # reversed because at the other end of the SMILES
    headgroup = headgroup_SMILES[lipidclass]
    #if lipidclass is a lyso PL and in sn2 position, modify headgroup as needed here
    # if len(sm_SMILES) == 2 and sn_SMILES[1] == "":
    # for lyso PL in sn2 position, modify headgroup as needed here
    COMPLETE_SMILES = sn1_SMILES + headgroup + sn2_SMILES[1:]

    return COMPLETE_SMILES


def get_lipid_SMILES(lipidclass, backbone, acyl_types, chains):
    sn_SMILES = list()
    #if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
        #ERROR

    for i in range(0, len(chains)):
        #print("C" * chains[i][0])
        if acyl_types[i] in ["acyl", "amide"]:
            sn_SMILES.append(add_dbs(("C" * chains[i][0]), acyl_db_table[chains[i][0]][chains[i][1]], "cis"))
        if acyl_types[i] in ["alkyl"]:
            sn_SMILES.append(add_dbs(("C" * chains[i][0]), alkyl_db_table[chains[i][0]][chains[i][1]], "cis"))
        if acyl_types[i] == "sbase":
            if chains[i][1] < 2: #db < 2
                sn_SMILES.append(add_dbs(("C" * chains[i][0]), sbase_db_table[chains[i][0]][chains[i][1]], "trans"))
            else:
                sn_SMILES.append(add_dbs(("C" * chains[i][0]), sbase_db_table[chains[i][0]][chains[i][1]][1:], "cis"))
                sn_SMILES.append(add_dbs(("C" * chains[i][0]), sbase_db_table[chains[i][0]][chains[i][1]][1], "trans"))

    if backbone == "bbSL":
        COMPLETE_SMILES = get_sphingolipid_SMILES(lipidclass, sn_SMILES)
    elif len(chains) < 3:
        COMPLETE_SMILES = get_glycerolipid_SMILES(lipidclass, sn_SMILES)

    elif:
        match lipidclass:
            case "TG":
                COMPLETE_SMILES = get_TG_SMILES(lipidclass, sn_SMILES)
            case "CL":
                COMPLETE_SMILES = get_CL_SMILES(lipidclass, sn_SMILES)
            case "LysoCL":
                COMPLETE_SMILES = get_LysoCL_SMILES(lipidclass, sn_SMILES)
            case "HemiBMP":
                COMPLETE_SMILES = get_HemiBMP_SMILES(lipidclass, sn_SMILES)
            case "BDP":
                COMPLETE_SMILES = get_BDP_SMILES(lipidclass, sn_SMILES)
    else:
        return "Not Available"

    #N-acyl_serine, N-acyl_ethanolamine, O-acylceramide all need to figure out SMILES.

    #print(COMPLETE_SMILES)
    return convert_to_canonical_SMILES(COMPLETE_SMILES)


#my_lipid = get_lipid_SMILES("Ceramide", "bbSL", ["sbase", "acyl"], [[18,1,1],[16,0,0]])
#print(my_lipid)
#make test set with pasted SMILES, then compare to the SMILES we get back from running this on the composition.  canonicalize both

