from rdkit import Chem
import copy

#this creates default positions for dbs.  specifics modified after
#first dimension in number of carbons, second is number of double bonds
acyl_db_table = {

}

for c in range(0,37):
    acyl_db_table[c] = {}
    for d in range(0,7):
        if (c > 5 and c < 22 and d > (c - 5) / 3): # or (c < 6 and d > 0):
            continue
        match d:  # #db
            case 0:
                acyl_db_table[c][d] = []
            case 1:
                if c < 11:
                    acyl_db_table[c][d] = [4]
                else:
                    acyl_db_table[c][d] = [9]
            case 2:
                if c < 16:
                    acyl_db_table[c][d] = [4,7]
                else:
                    acyl_db_table[c][d] = [9 ,12]
            case 3:
                if c < 18:
                    acyl_db_table[c][d] = [4, 7, 10]
                else:
                    acyl_db_table[c][d] = [9 ,12, 15]
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
for c in range(14,37):
    for d in range(1,7):
        if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
            continue
        if d == 1:
            alkyl_db_table[c][d] = [1]
        elif d > 1:
            alkyl_db_table[c][d] =[1] + acyl_db_table[c][d-1]


sbase_db_table = copy.deepcopy(acyl_db_table)
for c in range(14, 37):
    for d in range(1, 7):
        if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
            continue
        if d == 1:
            sbase_db_table[c][d] = [1]
        elif d > 1:
            sbase_db_table[c][d] = [1] + [x - 3 for x in acyl_db_table[c][d-1]]


def convert_to_random_smiles(smiles_str):
    random_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_str), isomericSmiles=True, doRandom=True)
    while(random_smiles[0] != 'C'  or random_smiles[1] != 'C' or random_smiles[-1] != 'C'):
        random_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_str), isomericSmiles=True, doRandom=True)
    return random_smiles


#take in any smiles and return canonical smiles
def convert_to_canonical_smiles(smiles_str):
    #print(smiles_str)
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

def add_hydroxyl(smile_str, position):
    smile_str_list = list(smile_str)
    if len(smile_str_list)  < position:
        return smile_str
    smile_str_list.insert(-(position-1), "(O)")
    new_smile_str = ''.join(smile_str_list)
    return new_smile_str

def reverse_acyl_smiles(smiles_str): #allows you to reverse the acyl groups so that they can point the other way, but maintain parentheses correctly
    new_smiles_str = ""
    for ch in smiles_str:
        if ch == "(":
            new_ch = ")"
        elif ch == ")":
            new_ch = "("
        else:
            new_ch = ch

        new_smiles_str = new_ch + new_smiles_str
    return new_smiles_str



def get_TG_smiles(lipidclass, sn_smiles):
    sn1_smiles = reverse_acyl_smiles(sn_smiles[0])[1:]
    sn2_smiles = reverse_acyl_smiles(sn_smiles[1])[1:]  # reversed because at the other end of the smiles
    sn3_smiles = reverse_acyl_smiles(sn_smiles[2])[1:]  # reversed because at the other end of the smiles
    complete_smiles = "O=C(OCC(OC(=O)" + sn1_smiles + ")COC(=O)" + sn2_smiles + ")" +  sn3_smiles
    return complete_smiles

def get_CL_smiles(lipidclass, sn_smiles):
    sn1_smiles = sn_smiles[0]
    sn2_smiles = reverse_acyl_smiles(sn_smiles[1])[1:]  # reversed because at the other end of the smiles
    sn3_smiles = reverse_acyl_smiles(sn_smiles[2])[1:]  # reversed because at the other end of the smiles
    sn4_smiles = reverse_acyl_smiles(sn_smiles[3])[1:]  # reversed because at the other end of the smiles
    complete_smiles = sn1_smiles + "(=O)OC[C@@H](OC(=O)" + sn2_smiles + ")COP(=O)(OCC(O)COP(=O)(OC[C@H](OC(=O)" + \
                        sn3_smiles + ")COC(=O)" + sn4_smiles + ")O)O"
    return complete_smiles


def get_LysoCL_smiles(lipidclass, sn_smiles):
    sn1_smiles = sn_smiles[0]
    sn2_smiles = reverse_acyl_smiles(sn_smiles[1])[1:]  # reversed because at the other end of the smiles
    sn3_smiles = reverse_acyl_smiles(sn_smiles[2])[1:]  # reversed because at the other end of the smiles

    complete_smiles = sn1_smiles + "(=O)OC[C@H](COP(=O)(O)OCC(COP(=O)(O)OC[C@@H](COC(=O)" + \
                        sn2_smiles + ")OC(=O)" + sn3_smiles + ")O)O"
    return complete_smiles


def get_BDP_smiles(lipidclass, sn_smiles):
    sn1_smiles = sn_smiles[0]
    sn2_smiles = reverse_acyl_smiles(sn_smiles[1])[1:]  # reversed because at the other end of the smiles
    sn3_smiles = reverse_acyl_smiles(sn_smiles[2])[1:]  # reversed because at the other end of the smiles
    sn4_smiles = reverse_acyl_smiles(sn_smiles[3])[1:]  # reversed because at the other end of the smiles
    complete_smiles = sn1_smiles + "(=O)OC[C@@H](COP(=O)(O)OC[C@H](COC(=O)" + sn2_smiles + ")OC(=O)" + \
                      sn3_smiles + ")OC(=O)" + sn4_smiles
    return complete_smiles


def get_HemiBMP_smiles(lipidclass, sn_smiles):
    sn1_smiles = sn_smiles[0]
    sn2_smiles = reverse_acyl_smiles(sn_smiles[1])[1:]  # reversed because at the other end of the smiles
    sn3_smiles = reverse_acyl_smiles(sn_smiles[2])[1:]  # reversed because at the other end of the smiles
    complete_smiles = sn1_smiles + "(=O)OC[C@H](O)COP(=O)(O)OC[C@@H](COC(=O)" + sn2_smiles + ")OC(=O)" + \
                      sn3_smiles
    return complete_smiles


def get_OAcylCeramide_smiles(lipidclass, sn_smiles):

    sn1_smiles = reverse_acyl_smiles(sn_smiles[0])[:-3] # reversed because at the other end of the smiles
    sn2_smiles = (sn_smiles[1])  # reversed because at the other end of the smiles
    sn3_smiles = reverse_acyl_smiles(sn_smiles[2])[:-1]

    complete_smiles = sn2_smiles + "(=O)OC[C@H](NC(=O)" + sn3_smiles + ")[C@H](O)" + \
                      sn1_smiles
    #print(complete_smiles)
    return complete_smiles



def get_N_Acyl_PE_smiles(lipidclass, sn_smiles):
    sn1_smiles = reverse_acyl_smiles(sn_smiles[0]) # reversed because at the other end of the smiles
    sn2_smiles = reverse_acyl_smiles(sn_smiles[1])  # reversed because at the other end of the smiles
    sn3_smiles = sn_smiles[2]

    complete_smiles = sn3_smiles + "(=O)NCCOP(O)(=O)OC[C@@H](COC(" + sn1_smiles[:-1] + ")=O)OC(=O)" + sn2_smiles[:-1]

    #"CCCCCCCCCCCCCCCC(=O)NCCOP(O)(=O)OC[C@@H](COC(CCCCC)=O)OC(=O)C"


    #complete_smiles = sn1_smiles + "(=O)OC[C@H](COP(O)(=O)OC[C@H](N" + sn3_smiles + ")C(O)=O)OC(=O)" + sn2_smiles[1:]

    return complete_smiles



def get_N_Acyl_PS_smiles(lipidclass, sn_smiles):
    sn1_smiles = reverse_acyl_smiles(sn_smiles[0]) # reversed because at the other end of the smiles
    sn2_smiles = reverse_acyl_smiles(sn_smiles[1])  # reversed because at the other end of the smiles
    sn3_smiles = sn_smiles[2]

    complete_smiles = sn3_smiles + "(=O)N[C@@H](CO)C(=O)COP(=O)(O)" + sn2_smiles + "(=O)OCC(O)" + \
                      sn1_smiles[:-1] + "(=O)CO"
    return complete_smiles



def get_sphingolipid_smiles(lipidclass, sn_smiles, chains):
    headgroup_smiles_1OH = { #for m LCBs, like m18:1
        # lipidclass:   #headgroup smiles

        "Ceramide": "[C@@H](O)[C@H](C)NC(=O)",
        "LCB": "[C@@H](O)[C@H](C)N",

    }

    headgroup_smiles_2OH = { #for d LCBs, like d18:1
        # lipidclass:   #headgroup smiles
        "AcGD1a": '[C@@H](O)[C@H](CO[C@H]1[C@@H]([C@@H](O)[C@@H]([C@H](O1)CO)O[C@@H]1O[C@H](CO)[C@H](O[C@@H]2O[C@H](CO)[C@@H]([C@H](O[C@H]3[C@@H]([C@H]([C@@H](O)[C@H](O3)CO)O[C@]3(OC([C@@H]([C@H](O)CO)O)[C@H](NC(C)=O)[C@@H](O)C3)C(O)=O)O)[C@H]2CC(=O)C)O)[C@@H]([C@H]1O)O[C@@]1(C[C@@H]([C@@H](NC(C)=O)C(O1)[C@H](O)[C@H](O)CO)O)C(O)=O)O)NC(=O)',
        "AcGD1b": '[C@@H](O)[C@H](CO[C@@H]1O[C@@H]([C@@H](O[C@H]2[C@H](O)[C@H]([C@H]([C@@H](CO)O2)O[C@H]2[C@H](CC(=O)C)[C@@H](O[C@H]3[C@@H]([C@H]([C@H]([C@H](O3)CO)O)O)O)[C@H]([C@@H](CO)O2)O)O[C@@]2(C[C@H](O)[C@@H](NC(C)=O)C([C@@H]([C@@H](CO)O[C@@]3(C[C@H](O)[C@H](C([C@H](O)[C@@H](CO)O)O3)NC(=O)C)C(=O)O)O)O2)C(O)=O)[C@H](O)[C@H]1O)CO)NC(=O)',
        "AcGD2": '[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3NC(C)=O)[C@H](O[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)C([C@H](O)[C@@H](CO)O[C@]4(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)C([C@H](O)[C@H](O)CO)O4)O3)[C@H]2O)[C@H](O)[C@H]1O)NC(=O)',
        "AcGD3": '[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)C([C@H](O)[C@@H](CO)O[C@]4(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)C([C@H](O)[C@H](O)CO)O4)O3)[C@H]2O)[C@H](O)C1O)NC(=O)',
        "AcGM1": '[C@@H](O)[C@H](CO[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O[C@@H]1O[C@H](CO)[C@@H]([C@H](O[C@]2(C(O)=O)C[C@H](O)[C@@H](NC(=O)C)[C@H]([C@@H]([C@H](O)CO)O)O2)[C@H]1O)O[C@H]1[C@@H]([C@H]([C@@H](O)[C@@H](CO)O1)O[C@H]1[C@H](O)[C@@H](O)[C@H]([C@H](O1)CO)O)NC(=O)C)O)O)NC(=O)',
        "AcGM2": '[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3NC(C)=O)[C@H](O[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)[C@H]([C@H](O)[C@H](O)CO)O3)[C@H]2O)[C@H](O)[C@H]1O)NC(=O)',
        "AcGM3": '[C@@H](O)[C@H](CO[C@@H]1OC(CO)[C@@H](O[C@@H]2OC(CO)[C@H](O)[C@H](O[C@]3(C(=O)O)CC(O)[C@@H](NC(C)=O)C([C@H](O)[C@H](O)CO)O3)C2O)[C@H](O)C1O)NC(=O)',
        "AcGQ1b": '[C@@H](O)[C@H](CO[C@H]1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O[C@H]3[C@H](NC(=O)C)[C@H]([C@@H](O)[C@H](O3)CO)O[C@@H]3O[C@H](CO)[C@@H]([C@H](O[C@@]4(C[C@H](O)[C@H]([C@@H](O4)[C@@H]([C@@H](CO)O[C@@]4(C[C@H](O)[C@H]([C@@H](O4)[C@H](O)[C@@H](CO)O)NC(C)=O)C(=O)O)O)NC(C)=O)C(O)=O)[C@H]3O)O)[C@H](O[C@]3(O[C@H]([C@@H]([C@@H](O)C3)NC(=O)C)[C@H](O)[C@@H](CO)O[C@@]3(C(=O)O)O[C@@H]([C@@H]([C@@H](CO)O)O)[C@@H]([C@@H](O)C3)NC(=O)C)C(=O)O)[C@H]2O)[C@H](O1)CO)O)NC(=O)',
        "AcGT1b": '[C@@H](O)[C@H](CO[C@H]1[C@H](O)[C@@H](O)[C@@H]([C@H](O1)CO)O[C@@H]1O[C@H](CO)[C@@H]([C@@H]([C@H]1O)O[C@]1(C(=O)O)C[C@H](O)[C@@H](NC(=O)C)C(O1)[C@H](O)[C@@H](CO)O[C@]1(C(O)=O)C[C@H](O)[C@@H](NC(C)=O)C(O1)[C@H](O)[C@H](O)CO)O[C@H]1[C@H](CC(=O)C)[C@H]([C@@H](O)[C@H](O1)CO)O[C@H]1[C@H](O)[C@@H](O[C@]2(C(O)=O)C[C@H](O)[C@@H](NC(=O)C)C(O2)[C@H](O)[C@H](O)CO)[C@@H](O)[C@@H](CO)O1)NC(=O)',
        "Ceramide": "[C@@H](O)[C@H](CO)NC(=O)",
        "Ceramide_P": "[C@@H](O)[C@H](COP(=O)(O)O)NC(=O)",
        "CPE": "[C@@H](O)[C@H](COP(=O)(O)OCCN)NC(=O)",
        "GB3": "C(O)C(COC1OC(CO)C(OC2OC(CO)C(OC3OC(CO)C(O)C(O)C3O)C(O)C2O)C(O)C1O)NC(=O)",
        "GcGM2": "[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3NC(C)=O)[C@H](O[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(CO)=O)[C@H]([C@H](O)[C@H](O)CO)O3)[C@H]2O)[C@H](O)[C@H]1O)NC(=O)", #thand Gemini
        "GcGM3": "[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(=O)CO)[C@H]([C@H](O)[C@H](O)CO)O3)[C@H]2O)[C@H](O)[C@H]1O)NC(=O)",
        "HexCer": "C(O)C(CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)",
        "LacCer": "C(O)C[C@H](O[C@H]1[C@@H](CO)[C@@H](O)[C@H](O)[C@@H](O[C@H]2[C@@H](CO)[C@@H](O)[C@H](O)[C@@H](O)O2)O1)NC(=O)",
        "LCB": "[C@@H](O)[C@@H](N)CO",
        "LCB_P": "[C@H]([C@H](COP(=O)(O)O)N)O",
        "LysoCPE": "[C@@H](O)[C@H](COP(=O)(O)OCCN)N",
        "LysoCPI":  "C(O)C(COP(O)(=O)OC1C(O)C(O)C(O)C(O)C1O)N",
        "LysoHexCer": "[C@H]([C@H](CO[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)N)O",
        "LysoLacCer": "[C@H]([C@H](CO[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O[C@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O)O)O)O)N)O",
        "LysoSM": "[C@H]([C@H](COP(=O)([O-])OCC[N+](C)(C)C)N)O",
        "SM": "[C@@H](O)[C@H](COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)",
        "Sulfatide": "[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)C(OS(=O)(=O)O)C1O)NC(=O)",
}

    headgroup_smiles_3OH = { #for t LCBs, like t18:0
        "Ceramide": "[C@@H](O)[C@@H](O)[C@H](CO)NC(=O)",
        "Ceramide_P": "(O)[C@@H](O)[C@H](COP(O)(O)=O)NC(=O)",
        "CPE": "C(O)[C@@H](O)[C@H](COP(=O)(O)OCCN)NC(=O)",
        "CPI": "C(O)[C@@H](O)[C@H](COP(OC1[C@@H]([C@@H](O)C([C@@H](O)[C@H]1O)O)O)(=O)O)NC(=O)",
        "LCB": "[C@H]([C@H]([C@H](CO)N)O)O",
        "LCB_P": "[C@H]([C@H]([C@H](COP(=O)(O)O)N)O)O",
        "LysoCPE": "C(O)[C@@H](O)[C@H](COP(=O)(O)OCCN)N",
        "LysoCPI": "C(O)[C@@H](O)[C@H](COP(OC1[C@@H]([C@@H](O)C([C@@H](O)[C@H]1O)O)O)(=O)O)N",
        "LysoHexCer": "C(O)[C@H]([C@H](CO[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)N)O",
        "MIP2C": "[C@@H](O)[C@@H](O)[C@H](COP(O)(=O)O[C@H]1[C@@H]([C@H]([C@H](O)[C@H]([C@H]1O)O[C@H]1[C@H]([C@@H](O)[C@@H]([C@@H](COP(=O)(O)OC2[C@@H]([C@@H](O)C([C@H]([C@H]2O)O)O)O)O1)O)O)O)O)NC(=O)",
        "MIPC": "C(O)[C@@H](O)[C@H](COP(=O)(O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)NC(=O)",
        "LysoSM": "C(O)[C@H]([C@H](COP(=O)([O-])OCC[N+](C)(C)C)N)O",

    }


#since LCB hydroxyls are incorporated in the head group portion of the SMILES, we need different headgroups for m,d,t LCBs
    if chains[0][2] == 0:
        headgroup_smiles = headgroup_smiles_1OH
    elif chains[0][2] == 1:
        headgroup_smiles = headgroup_smiles_2OH
    elif chains[0][2] == 2:
        headgroup_smiles = headgroup_smiles_3OH
    else:
        return("")

    if lipidclass not in headgroup_smiles:
        return ""


    sn1_smiles = sn_smiles[0]
    if (chains[0][2] == 2):
        complete_smiles = sn1_smiles[4:] + headgroup_smiles[lipidclass]
    else:
        complete_smiles = sn1_smiles[3:] + headgroup_smiles[lipidclass]

    if (len(sn_smiles) == 2):
        sn2_smiles = reverse_acyl_smiles(sn_smiles[1])  # reversed because at the other end of the smiles
        complete_smiles += sn2_smiles[1:]

    return complete_smiles

#currently does not work for Lysolipids with acyl at SN2 position
def get_glycerolipid_smiles(lipidclass, sn_smiles):

    # Headgroups here are different from that for fragmentation.
    # They generally go from the first non-carbon atom to the last non-carbon atom, or parenthesis/bracket
    headgroup_smiles = {
        # lipidclass:   #headgroup smiles
        #glycero and phospholipids
        "Alkyl_LPA": "[C@](COP(=O)(O)O)([H])(O)CO",
        "Alkyl_LPC": "OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)O",
        "Alkyl_LPE": "OC[C@H](COP(=O)(O)OCCN)O",
        "Alkyl_LPG": "[C@]([H])(O)(COP(=O)(O)OC[C@@]([H])(O)CO)CO",
        "Alkyl_LPS": "OC[C@H](COP(=O)(O)OC[C@@H](C(=O)O)N)O",
        "Alkyl_PE": "OC[C@H](COP(=O)(O)OCCN)OC(=O)",
        "Alkyl_PC": "OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)",
        "Alkyl_PS": "OC[C@H](COP(=O)(O)OC[C@@H](C(=O)O)N)OC(=O)",
        "BMP": "(=O)OC[C@@H](O)COP(=O)(O)OC[C@@H](O)COC(=O)",
        "CDP_DG": "(=O)OC[C@H](COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@@H](C([C@@H](O1)N2C=CC(=NC2=O)N)O)O)OC(=O)",
        "DG": "(=O)OC[C@H](CO)OC(=O)",
        "DGDG": "(=O)OC[C@H](CO[C@@H]1O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)C(O)C2O)[C@H](O)C(O)C1O)OC(=O)",
        "DGTS": "(=O)OC[C@H](COCCC(C(=O)[O-])[N+](C)(C)C)OC(=O)",
        "DMPE": "(=O)OC[C@H](COP(=O)(O)OCCN(C)C)OC(=O)",
        "FA": "(=O)O",
        "LPA": "(=O)OCC(COP(=O)(O)O)O",
        "LPC": "(=O)OC[C@@H](COP(=O)([O-])OCC[N+](C)(C)C)O",
        "LPE": "(=O)OC[C@H](COP(=O)(O)OCCN)O",
        "LPG": "(=O)OC[C@@H](COP(=O)(O)OC[C@H](CO)O)O",
        "LPI": "(=O)OC[C@H](COP(=O)(O)OC1[C@@H]([C@H](C([C@H]([C@H]1O)O)O)O)O)O",
        "LPS": "(=O)OC[C@H](COP(=O)(O)OC[C@@H](C(=O)O)N)O",
        "MG": "(=O)OCC(CO)O",
        "MGDG": "(=O)OC[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)OC(=O)",
        "MMPE": "(=O)OC[C@H](COP(=O)(O)OCCNC)OC(=O)",
        "PA": "(=O)OCC(COP(=O)(O)O)OC(=O)",
        "PC": "(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)",
        "PE": "(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)",
        "PG": "(=O)OC[C@H](COP(=O)(O)OC[C@H](CO)O)OC(=O)",
        "PI": "(=O)OC[C@H](COP(=O)(O)OC1[C@@H]([C@H](C([C@H]([C@H]1O)O)O)O)O)OC(=O)",
        "PS": "(=O)OC[C@H](COP(=O)(O)OC[C@@H](C(=O)O)N)OC(=O)",

        #the rest
        "APCS": "(=O)N[C@H](C(=O)O)COP(=O)([O-])OCC[N+](C)(C)C",
        "Carn": "(=O)O[C@H](CC(=O)[O-])C[N+](C)(C)C",
        "CE": "(=O)O[C@H]1CC[C@@]2([C@H]3CC[C@]4([C@H]([C@@H]3CC=C2C1)CC[C@@H]4[C@H](C)CCCC(C)C)C)C",
        "ErgE": "(=O)O[C@H]1CC[C@@]2([C@H]3CC[C@]4([C@H](C3=CC=C2C1)CC[C@@H]4[C@H](C)/C=C/[C@H](C)C(C)C)C)C",
        "Ethanolamine": "(=O)NCCO",
        "FA": "(=O)O",
        "Taurine": "(=O)NCCS(=O)(=O)O",
    }

    sn2_lyso_headgroup_smiles = {
        # lipidclass:   #headgroup smiles
    "LPA": "(=O)OCC(COP(=O)(O)O)O",
    "LPA": "OCC(COP(=O)(O)O)OC(=O)",
    "LPC": "OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)",
    "LPE": "OC[C@H](COP(=O)(O)OCCN)OC(=O)",
    "LPG": "OC[C@@H](COP(=O)(O)OC[C@H](CO)O)OC(=O)",
    "LPI": "OC[C@H](COP(=O)(O)OC1[C@@H]([C@H](C([C@H]([C@H]1O)O)O)O)O)OC(=O)",
    "LPS": "OC[C@H](COP(=O)(O)OC[C@@H](C(=O)O)N)OC(=O)"
    }

    if lipidclass not in headgroup_smiles:
        return ""

    # if lipidclass is a lyso PL and in sn2 position,
    if(len(sn_smiles) == 2 and sn_smiles[0] == ""):
        headgroup = sn2_lyso_headgroup_smiles[lipidclass]
    else: # all others
        headgroup = headgroup_smiles[lipidclass]

    sn1_smiles = sn_smiles[0]
    if (len(sn_smiles) > 1):
        sn2_smiles = sn_smiles[1][::-1]  # reversed because at the other end of the smiles
        complete_smiles = sn1_smiles + headgroup + sn2_smiles[1:]
    else:
        complete_smiles = sn1_smiles + headgroup





    # if len(sm_smiles) == 2 and sn_smiles[1] == "":
    # for lyso PL in sn2 position, modify headgroup as needed here

    return complete_smiles


def get_lipid_smiles(lipidclass, backbone, acyl_types, chains):
    sn_smiles = list()
    #print(lipidclass, chains)
    for i in range(0, len(chains)):
        if (chains[i][0] > 5 and chains[i][0] < 22 and chains[i][1] > (chains[i][0] - 5) / 3) or (chains[i][0] < 6 and chains[i][1] > 0):
            return ("")



        if acyl_types[i] in ["acyl", "amide", "alkyl"] and chains[i][2] > 1:
            return ""
        if acyl_types[i] in ["acyl", "amide"]:
            new_smiles = (add_dbs(("C" * chains[i][0]), acyl_db_table[chains[i][0]][chains[i][1]], "cis"))
            if chains[i][2] == 1:
               # if (chains[i][0] < 4):
               #     return ""
                if (i ==0) :
                    sn_smiles.append(add_hydroxyl(new_smiles, 2)) #assuming hydroxyl on acyl is at position 2
                else:
                    sn_smiles.append(add_hydroxyl(new_smiles, 3)) #becasue one c gets chopped later

            else:
                sn_smiles.append(new_smiles)
        if acyl_types[i] in ["alkyl"]:
            sn_smiles.append(add_dbs(("C" * chains[i][0]), alkyl_db_table[chains[i][0]][chains[i][1]], "cis"))

        if acyl_types[i] == "sbase":
            if chains[i][1] < 2: #db < 2
                sn_smiles.append(add_dbs(("C" * chains[i][0]), sbase_db_table[chains[i][0]][chains[i][1]], "trans"))
            else:
                cis_bonds = add_dbs(("C" * chains[i][0]), sbase_db_table[chains[i][0]][chains[i][1]][1:], "cis")
                sn_smiles.append(add_dbs((cis_bonds), sbase_db_table[chains[i][0]][chains[i][1]][0:1], "trans"))
    if backbone == "bbSL" and len(chains) < 3:
        complete_smiles = get_sphingolipid_smiles(lipidclass, sn_smiles, chains)
    elif len(chains) < 3:
        complete_smiles = get_glycerolipid_smiles(lipidclass, sn_smiles)

    else: #3 or more chains
        match lipidclass:
            case "TG":
                complete_smiles = get_TG_smiles(lipidclass, sn_smiles)
            case "CL":
                complete_smiles = get_CL_smiles(lipidclass, sn_smiles)
            case "LysoCL":
                complete_smiles = get_LysoCL_smiles(lipidclass, sn_smiles)
            case "HemiBMP":
                complete_smiles = get_HemiBMP_smiles(lipidclass, sn_smiles)
            case "BDP":
                complete_smiles = get_BDP_smiles(lipidclass, sn_smiles)
            case "OAcylCeramide":
               complete_smiles = get_OAcylCeramide_smiles(lipidclass, sn_smiles)
            case "N_Acyl_PE":
                 complete_smiles = get_N_Acyl_PE_smiles(lipidclass, sn_smiles)
            case "N_Acyl_PS":
                 complete_smiles = get_N_Acyl_PS_smiles(lipidclass, sn_smiles)
            case _:
                return ""

    return convert_to_canonical_smiles(complete_smiles)

#for testing:
#my_lipid = get_lipid_smiles("OAcylCeramide", "bbAcylSL", ["sbase", "amide", "acyl"], [[14,0,1],[2,0,0],[2,0,0]])
#print(my_lipid)

