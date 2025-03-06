import re
import sys


H2O = 18.01056
NH3 = 17.00274
PROTON = 1.00727645229
HYDROGEN = 1.0078250321
OXYGEN = 15.9949146221
ELECTRON = 0.0005489

ADDUCT = {
    "[M]": 0.0,
    "[M+H]+": 1.00727645229,
    "[M+Li]+": 7.016003,
    "[M+Na]+": 22.989218,
    "[M+NH4]+": 18.033823,
    "[M+K]+": 38.963708,
    "[M-H]-": -1.00727645229,
    "[M+Cl]-": 34.9694016,
    "[M+37Cl]-": 36.9664516,
    "[M+FA-H]-": 44.99820284,
    "[M-CH3]-": -15.0234750963,
    "[M+AcOH-H]-": 59.0138529,
    "[M+2H]2+": 2.0145529046,
    "[M-2H]2-": -2.0145529046,
    "[M-3H]3-": -3.0218293568,
}

ADDUCT_CHARGE = {
    "[M]": 1,
    "[M+H]+": 1,
    "[M+Li]+": 1,
    "[M+Na]+": 1,
    "[M+NH4]+": 1,
    "[M+K]+": 1,
    "[M-H]-": -1,
    "[M+Cl]-": -1,
    "[M+37Cl]-": -1,
    "[M+FA-H]-": -1,
    "[M-CH3]-": -1,
    "[M+AcOH-H]-": -1,
    "[M+2H]2+": 2,
    "[M-2H]2-": -2,
    "[M-3H]3-": -3,
}


POS_ADDUCTS = {"[M]", "[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+", "[M+2H]2+", "[M+Li]+"}
NEG_ADDUCTS = {
    "[M-H]-",
    "[M+Cl]-",
    "[M+37Cl]-",
    "[M+FA-H]-",
    "[M-CH3]-",
    "[M+AcOH-H]-",
    "[M-2H]2-",
    "[M-3H]3-",
}

ALL_ADDUCTS = POS_ADDUCTS.copy()
ALL_ADDUCTS |= NEG_ADDUCTS.copy()

ALL_LIPID_CLASSES = {
    "AcGD1a",
    "AcGD1b",
    "AcGD2",
    "AcGD3",
    "AcGM1",  # 5
    "AcGM2",
    "AcGM3",
    "AcGQ1b",
    "AcGT1b",
    "Alkyl_LPC",  # 10
    "Alkyl_LPE",
    "Alkyl_LPS",
    "Alkyl_PC",
    "Alkyl_PE",
    "Alkyl_PS",  # 15
    "APCS",
    "BDP",
    "BMP",
    "Carn",
    "CDP_DG",  # 20
    "CE",
    "Ceramide",
    "Ceramide_P",
    "CL",
    "CPE",  # 25
    "CPI",
    "DG",
    "DGDG",
    "DGTS",
    "DMPE",  # 30
    "ErgE",
    "Ethanolamine",
    "FA",
    "FAHFA",
    "GB3",  # 35
    "GcGM2",
    "GcGM3",
    "HemiBMP",
    "HexCer",
    "LacCer",  # 40
    "LCB_P",
    "LCB",
    "LPA",
    "LPC",
    "LPE",  # 45
    "LPG",
    "LPI",
    "LPS",
    "LysoCL",
    "LysoCPE",  # 50
    "LysoCPI",
    "LysoHexCer",
    "LysoSM",
    "MG",
    "MGDG",  # 55
    "MIP2C",
    "MIPC",
    "MMPE",
    "N_Acyl_PE",
    "N_Acyl_PS",  # 60
    "PA",
    "PC",
    "PE",
    "PG",
    "PI",  # 65
    "PS",
    "SM",
    "Sulfatide",
    "Taurine",
    "TG",  # 70
    "OAcylCeramide",
    "RetinylE"
}


def MW(form):  # calc MW from string that is chemical formula
    parsed = re.findall(r"([A-Z][a-z]*)(\d*)", form)
    weight = 0.00

    parsed = [list(elem) for elem in parsed]  # convert tuples from re to lists
    for i in range(0, len(parsed)):
        if parsed[i][1] == "":
            parsed[i][1] = 1  # set elements with no number, because only one, to one
    for elem in parsed:
        if elem[0] == "C":
            weight += 12 * int(elem[1])
        elif elem[0] == "H":
            weight += 1.0078250321 * int(elem[1])
        elif elem[0] == "N":
            weight += 14.0030740052 * int(elem[1])
        elif elem[0] == "O":
            weight += 15.9949146221 * int(elem[1])
        elif elem[0] == "P":
            weight += 30.97376151 * int(elem[1])
        elif elem[0] == "S":
            weight += 31.972070 * int(elem[1])

    return weight


def MW_list(atoms=list()):  # calc MW from list that is C,H,N,O,P,S
    (nC, nH, nN, nO, nP, nS) = atoms
    return (
        12.000 * nC
        + 1.0078250321 * nH
        + 14.0030740052 * nN
        + 15.9949146221 * nO
        + 30.97376151 * nP
        + 31.972070 * nS
    )


# only works on fatty acyls
def NL(chain):
    nC = chain[0]
    nH = chain[0] * 2 - (2 * chain[1])
    nO = 2 + chain[2]
    return MW_list([nC, nH, 0, nO, 0, 0])


# def LCB(nC,nD,hydroxyls):
#     nH = nC*2+3-(2*nD)
#     nO = 1 + hydroxyls
#     nN = 1
#     return MW_list([nC,nH,nN,nO,0,0])


LIPID_BACKBONES = {
    # backbone,   #linkage
    "PC": ["bbGPL", ["acyl", "acyl"]],
    "LPC": ["bbGPL", ["acyl", "acyl"]],
    "LPS": ["bbGPL", ["acyl", "acyl"]],
    "LPG": ["bbGPL", ["acyl", "acyl"]],
    "LPA": ["bbGPL", ["acyl", "acyl"]],
    "LPI": ["bbGPL", ["acyl", "acyl"]],
    "LPE": ["bbGPL", ["acyl", "acyl"]],
    "PE": ["bbGPL", ["acyl", "acyl"]],
    "MMPE": ["bbGPL", ["acyl", "acyl"]],
    "DMPE": ["bbGPL", ["acyl", "acyl"]],
    "PS": ["bbGPL", ["acyl", "acyl"]],
    "PI": ["bbGPL", ["acyl", "acyl"]],
    "PA": ["bbGPL", ["acyl", "acyl"]],
    "PG": ["bbGPL", ["acyl", "acyl"]],
    "SM": ["bbSL", ["sbase", "amide"]],
    "LysoSM": ["bbSL", ["sbase"]],
    "APCS": ["bbCarn", ["amide"]],
    "LCB": ["bbSL", ["sbase"]],
    "Alkyl_PS": ["bbGPL", ["alkyl", "acyl"]],
    "Alkyl_PE": ["bbGPL", ["alkyl", "acyl"]],
    "Alkyl_PC": ["bbGPL", ["alkyl", "acyl"]],
    "Alkyl_LPC": ["bbGPL", ["alkyl"]],
    "Alkyl_LPE": ["bbGPL", ["alkyl"]],
    "Alkyl_LPS": ["bbGPL", ["alkyl"]],
    "Ceramide": ["bbSL", ["sbase", "amide"]],
    "HexCer": ["bbSL", ["sbase", "amide"]],
    "LacCer": ["bbSL", ["sbase", "amide"]],
    "GB3": ["bbSL", ["sbase", "amide"]],
    "LysoHexCer": ["bbSL", ["sbase"]],
    "Sulfatide": ["bbSL", ["sbase", "amide"]],
    "Ceramide_P": ["bbSL", ["sbase", "amide"]],
    "MGDG": ["bbGL", ["acyl", "acyl"]],
    "DGDG": ["bbGL", ["acyl", "acyl"]],
    "TG": ["bbGL", ["acyl", "acyl", "acyl"]],
    "DG": ["bbGL", ["acyl", "acyl"]],
    "MG": ["bbGL", ["acyl"]],
    "CE": ["bbSterol", ["acyl"]],
    "ErgE": ["bbSterol", ["acyl"]],
    "RetinylE": ["bbRetinol", ["acyl"]],
    "Carn": ["bbCarn", ["acyl"]],
    "FA": ["bbFatty", ["acyl"]],
    "CL": ["bbCL", ["acyl", "acyl", "acyl", "acyl"]],
    "LysoCL": ["bbCL", ["acyl", "acyl", "acyl"]],
    "AcGM3": ["bbSL", ["sbase", "amide"]],
    "AcGM2": ["bbSL", ["sbase", "amide"]],
    "AcGM1": ["bbSL", ["sbase", "amide"]],
    "AcGD1a": ["bbSL", ["sbase", "amide"]],
    "AcGD1b": ["bbSL", ["sbase", "amide"]],
    "AcGT1b": ["bbSL", ["sbase", "amide"]],
    "AcGQ1b": ["bbSL", ["sbase", "amide"]],
    "AcGD2": ["bbSL", ["sbase", "amide"]],
    "AcGD3": ["bbSL", ["sbase", "amide"]],
    "GcGM3": ["bbSL", ["sbase", "amide"]],
    "GcGM2": ["bbSL", ["sbase", "amide"]],
    "N_Acyl_PS": ["bbAcylGPL", ["acyl", "acyl", "amide"]],
    "N_Acyl_PE": ["bbAcylGPL", ["acyl", "acyl", "amide"]],
    "Ethanolamine": ["bbCarn", ["amide"]],
    "CPI": ["bbSL", ["sbase", "amide"]],
    "MIPC": ["bbSL", ["sbase", "amide"]],
    "MIP2C": ["bbSL", ["sbase", "amide"]],
    "CoQ": ["bbIsoprenoid", ["isoprenoid"]],
    "BMP": ["bbGPL", ["acyl", "acyl"]],
    "Taurine": ["bbCarn", ["amide"]],
    "FAHFA": ["bbFAHFA", ["acyl", "acyl"]],
    "CPE": ["bbSL", ["sbase", "amide"]],
    "N_methylLCB": ["bbSL", ["sbase"]],
    "N_dimethylLCB": ["bbSL", ["sbase"]],
    "DGTS": ["bbGL", ["acyl", "acyl"]],
    "LCB_P": ["bbSL", ["sbase"]],
    "LysoCPI": ["bbSL", ["sbase"]],
    "LysoCPE": ["bbSL", ["sbase"]],
    "HemiBMP": ["bbGPL", ["acyl", "acyl", "acyl"]],
    "BDP": ["bbGPL", ["acyl", "acyl", "acyl", "acyl"]],
    "CDP_DG": ["bbGPL", ["acyl", "acyl"]],
    "OAcylCeramide": ["bbAcylSL", ["sbase", "amide", "acyl"]],

}

LIPID_HEADS = {
    # C     #H       #N    #O    #P    #S
    "PC": [5 + 3, 15 + 3, 1, 4, 1, 0],
    "LPC": [5 + 3, 15 + 3, 1, 4, 1, 0],
    "Alkyl_PE": [5, 12, 1, 4, 1, 0],
    "Alkyl_LPE": [5, 12, 1, 4, 1, 0],
    "PE": [5, 12, 1, 4, 1, 0],
    "MMPE": [6, 14, 1, 4, 1, 0],
    "DMPE": [7, 16, 1, 4, 1, 0],
    "PS": [6, 12, 1, 6, 1, 0],
    "Alkyl_PS": [6, 12, 1, 6, 1, 0],
    "Alkyl_LPS": [6, 12, 1, 6, 1, 0],
    "LPS": [6, 12, 1, 6, 1, 0],
    "LPG": [6, 13, 0, 6, 1, 0],
    "PI": [6 + 3, 12 + 5, 0, 9, 1, 0],
    "LPI": [6 + 3, 12 + 5, 0, 9, 1, 0],
    "PA": [3, 2 + 5, 0, 4, 1, 0],
    "PG": [6, 13, 0, 6, 1, 0],
    "LPA": [3, 2 + 5, 0, 4, 1, 0],
    "LPE": [5, 12, 1, 4, 1, 0],
    "SM": [5, 13, 1, 4, 1, 0],
    "LysoSM": [5, 13, 1, 4, 1, 0],
    "APCS": [8, 18, 1, 6, 1, 0],
    "LCB": [0, 1, 0, 1, 0, 0],
    "Alkyl_PC": [5 + 3, 15 + 3, 1, 4, 1, 0],
    "Alkyl_LPC": [5 + 3, 15 + 3, 1, 4, 1, 0],
    "Ceramide": [0, 1, 0, 1, 0, 0],
    "HexCer": [6, 11, 0, 6, 0, 0],
    "LacCer": [12, 21, 0, 11, 0, 0],
    "GB3": [18, 31, 0, 16, 0, 0],
    "LysoHexCer": [6, 11, 0, 6, 0, 0],
    "Sulfatide": [6, 11, 0, 9, 0, 1],
    "Ceramide_P": [0, 2, 0, 4, 1, 0],
    "MGDG": [6 + 3, 10 + 5, 0, 5, 0, 0],
    "DGDG": [12 + 3, 20 + 5, 0, 10, 0, 0],
    "TG": [3, 5, 0, 0, 0, 0],
    "DG": [3, 5, 0, 0, 0, 0],
    "MG": [3, 5, 0, 0, 0, 0],
    "CE": [27, 45, 0, 0, 0, 0],
    "RetinylE": [20,29,0,0,0,0],
    "ErgE": [28, 43, 0, 0, 0, 0],
    "Carn": [7, 14, 1, 2, 0, 0],
    "FA": [0, 1, 0, 0, 0, 0],
    "CL": [9, 18, 0, 9, 2, 0],
    "LysoCL": [9, 18, 0, 9, 2, 0],
    "AcGM3": [23, 38, 1, 19, 0, 0],
    "AcGM2": [23 + 8, 38 + 15 - 2, 1 + 1, 19 + 6 - 1, 0, 0],
    "AcGM1": [23 + 8 + 6, 38 + 15 - 2 + 10, 1 + 1, 19 + 6 + 5 - 1, 0, 0],
    "AcGD2": [48 - 6, 78 - 10, 3, 37 - 5, 0, 0],
    "AcGD3": [48 - 6 - 8, 78 - 10 - 13, 3 - 1, 37 - 5 - 5, 0, 0],
    "AcGD1a": [48, 78, 3, 37, 0, 0],
    "AcGD1b": [48, 78, 3, 37, 0, 0],
    "AcGT1b": [48 + 11, 78 + 17, 3 + 1, 37 + 8, 0, 0],
    "AcGQ1b": [48 + 22, 78 + 34, 3 + 2, 37 + 16, 0, 0],
    "GcGM3": [23, 38, 1, 20, 0, 0],
    "GcGM2": [23 + 8, 38 + 15 - 2, 1 + 1, 20 + 6 - 1, 0, 0],
    "N_Acyl_PS": [6, 12 - 1, 0, 6, 1, 0],
    "N_Acyl_PE": [5, 12 - 1, 0, 4, 1, 0],
    "Ethanolamine": [2, 5, 0, 1, 0, 0],
    "CPI": [6, 12, 0, 9, 1, 0],
    "MIPC": [12, 22, 0, 14, 1, 0],
    "MIP2C": [18, 33, 0, 22, 2, 0],
    "CoQ": [9, 9, 0, 4, 0, 0],
    "BMP": [6, 13, 0, 6, 1, 0],
    "Taurine": [2, 6, 0, 3, 0, 1],
    "FAHFA": [0, 1, 0, 0, 0, 0],
    "CPE": [2, 7, 1, 4, 1, 0],
    "N_methylLCB": [1, 3, 0, 1, 0, 0],
    "N_dimethylLCB": [2, 5, 0, 1, 0, 0],
    "DGTS": [10, 18, 1, 2, 0, 0],
    "LCB_P": [0, 2, 0, 4, 1, 0],
    "LysoCPI": [6, 12, 0, 9, 1, 0],
    "LysoCPE": [2, 7, 1, 4, 1, 0],
    "HemiBMP": [6, 12, 0, 5, 1, 0],
    "BDP": [6, 11, 0, 4, 1, 0],
    "CDP_DG": [12, 19, 3, 11, 2, 0],
    "OAcylCeramide": [0,0,0,0,0,0]
}


LIPID_MAPS = {
    "PC": "GP0101",
    "LPC": "GP0105",
    "Alkyl_PE": "GP0202",
    "Alkyl_LPE": "GP0206/GP0207",
    "PE": "GP0201",
    "MMPE": "GP0201",
    "DMPE": "GP0201",
    "PS": "GP0301",
    "Alkyl_PS": "GP0202",
    "Alkyl_LPS": "GP0206/GP0207",
    "LPS": "GP0205",
    "LPG": "GP0405",
    "PI": "GP0601",
    "LPI": "GP0605",
    "PA": "GP1001",
    "PG": "GP0401",
    "LPA": "GP1005",
    "LPE": "GP0205",
    "SM": "SP0301",
    "APCS": "???",
    "LysoSM": "SP0106",
    "LCB": "SP0101/SP0102/SP0103",
    "Alkyl_PC": "GP0102",
    "Alkyl_LPC": "GP0106/GP0107",
    "Ceramide": "SP201/SP202/SP203",
    "HexCer": "SP0501",
    "LacCer": "SP0501",
    "GB3": "SP0502 ",
    "LysoHexCer": "SP0106",
    "Sulfatide": "SP0602",
    "Ceramide_P": "SP0205",
    "MGDG": "GL0501",
    "DGDG": "GL0501",
    "TG": "GL0301",
    "DG": "GL0201",
    "MG": "GL0101",
    "CE": "ST0101",
    "ErgE": "??",
    "Carn": "FA0707",
    "FA": "FA0101",
    "CL": "GP1201",
    "LysoCL": "GP1202",
    "AcGM3": "SP0601",
    "AcGM2": "SP0601",
    "AcGM1": "SP0601",
    "AcGD2": "SP0601",
    "AcGD3": "SP0601",
    "AcGD1a": "SP0601",
    "AcGD1b": "SP0601",
    "AcGT1b": "SP0601",
    "AcGQ1b": "SP0601",
    "GcGM3": "SP0601",
    "GcGM2": "SP0601",
    "N_Acyl_PS": "??",
    "N_Acyl_PE": "??",
    "Ethanolamine": "FA0804",
    "CPI": "SP0303",
    "MIPC": "SP0303",
    "MIP2C": "SP0303",
    "CoQ": "PR0201",
    "BMP": "GP0410",
    "Taurine": "FA0804",
    "FAHFA": "FA0701",
    "CPE": "SP0302",
    "N_methylLCB": "SP0107",
    "N_dimethylLCB": "SP0107",
    "DGTS": "GL00",
    "LCB_P": "SP0105",
    "LysoCPI": "SP0303",
    "LysoCPE": "SP0302",
    "HemiBMP": "GP0409",
    "BDP": "GP0408",
    "CDP_DG": "GP1301",
    "OAcylCeramide": "SP0204",
    "RetinylE": "PR0109"
}


class Lipid:
    def __init__(self):
        self.lipidclass = ""
        self.chains = list()
        self.adduct = ""
        self.lipidName = ""


    # complete set of adducts, for future reference for DI libraries
    #neg_adduct_set = ["[M-H]-", "[M+Cl]-", "[M+37Cl]-", "[M+FA-H]-", "[M+AcOH-H]-"]
    #pos_adduct_set = ["[M+Li]+", "[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+"]

    #current set of adducts we wish to search for
    neg_adduct_set = ["[M-H]-","[M+FA-H]-"]
    pos_adduct_set = ["[M+H]+", "[M+Na]+", "[M+NH4]+"]
    NCE = "20,30,40"

    # lipidclass is one of the above listed lipid backbones
    # chains is a list of lists, where each member list is [#carbons, #double bond, #hydroxyls] for one carbon chain in the lipid
    def set_chains_and_adduct(self, lipidclass, chains=[], adduct="", lipidName=""):
        self.lipidclass = lipidclass
        self.chains = chains
        self.adduct = adduct
        self.lipidName = lipidName

        if not self.adduct:
            raise Exception("Please define adduct")

        if self.adduct not in ADDUCT:
            raise Exception("Unknown adduct" + self.adduct)

        if lipidclass not in LIPID_BACKBONES:
            raise Exception("Missing lipid class" + lipidName)

        if lipidclass not in LIPID_MAPS:
            raise Exception("Missing lipid maps classifier for" + lipidName)

    def MF(self):
        formula = list(LIPID_HEADS[self.lipidclass])  # clone head group
        bbType = LIPID_BACKBONES[self.lipidclass][0]
        chainType = LIPID_BACKBONES[self.lipidclass][1]
        C = 0
        H = 1
        N = 2
        O = 3
        P = 4
        S = 5
        if bbType != "bbIsoprenoid":
            for i in range(0, len(chainType)):
                chain = self.chains[i]
                formula[C] += chain[0]

                if chainType[i] == "acyl":
                    formula[O] += 2
                    formula[H] += (chain[0] * 2) - 1

                elif chainType[i] == "alkyl":
                    formula[O] += 1
                    formula[H] += (chain[0] * 2) + 1

                elif chainType[i] == "amide":
                    formula[H] += (chain[0] * 2) - 1
                    formula[N] += 1
                    formula[O] += 1  # newly added for bug

                elif chainType[i] == "sbase":
                    formula[H] += (chain[0] * 2) + 1

                else:
                    formula[H] += (chain[0] * 2) - 1

                formula[H] -= chain[1] * 2  # subtract two H for each double bond
                formula[O] += chain[2]  # add an O for each hydroxyl group

            if len(chainType) == 1:  # for Lyso or MG
                if bbType == "bbGPL":
                    formula[O] += 1
                    formula[H] += 1
                elif bbType == "bbSL":
                    formula[N] += 1
                    formula[H] += 1
                elif bbType == "bbGL":
                    formula[O] += 2
                    formula[H] += 2

            if (
                len(chainType) == 2
                and bbType == "bbGPL"
                and (self.chains[0][0] == 0 or self.chains[1][0] == 0)
            ):  # for Lyso PLs
                formula[O] -= 1
                formula[H] += 2

            if len(chainType) == 2 and bbType == "bbGL":  # for DG
                formula[O] += 1
                formula[H] += 1



            if len(chainType) < 4 and bbType == "bbCL":
                lysoCount = 4 - len(chainType)
                formula[O] += lysoCount
                formula[H] += lysoCount

            if bbType == "bbFAHFA":
                formula[H] -= 1  # loss of hydrogen from where the second FA is attached

        else:  # for isoprenoids only, where each # = isoprenyl units
            formula[C] += self.chains[0][0] * 5
            formula[H] += self.chains[0][0] * 8 + 1  # +! is for terminal isoprene unit

        return formula

    def prettyFormula(self):
        atoms = self.MF()
        atomNames = ["C", "H", "N", "O", "P", "S"]

        formulaStr = ""
        for i in range(0, len(atoms)):
            if atoms[i] >= 1:
                formulaStr += atomNames[i]
            if atoms[i] > 1:
                formulaStr += str(atoms[i])

        return formulaStr

    def printNist(self):
        bbType = LIPID_BACKBONES[self.lipidclass][0]
        Nchains = len(self.chains)
        totalC = 0
        totalD = 0
        h = []
        for chain in self.chains:
            h.append(chain[2])  # list of hydroxyls on each carbon chain
        totalH = sum(h)
        ChainStrings = []
        if bbType != "bbIsoprenoid":
            for i in range(0, Nchains):
                totalC += self.chains[i][0]
                totalD += self.chains[i][1]

                d = self.chains[i][1]
                chainstr = ""
                if (
                    LIPID_BACKBONES[self.lipidclass][1][i] == "alkyl"
                ):  # to name plasmalogens ets with p and o convention
                    if d > 0:
                        chainstr += "p-"
                        d -= 1
                    else:
                        chainstr += "o-"
                chainstr += str(self.chains[i][0])  # of carbons
                if len(self.chains) > i:
                    chainstr += ":" + str(d)  # double bonds
                if (
                    i == 0 and LIPID_BACKBONES[self.lipidclass][1][0] == "sbase"
                ):  # if sn1 is a sphingoid base
                    h[0] += 1
                if h[i] > 0:
                    chainstr += ";O"
                if h[i] > 1:
                    chainstr += str(h[i])
                ChainStrings.append(chainstr)

            if LIPID_BACKBONES[self.lipidclass][1][0] == "sbase":
                totalH += 1
            SumName = (
                re.sub("Alkyl_", "", self.lipidclass) + "("
            )  # to removed Alkyl labels from plasmalogens
            SumName += str(totalC) + ":" + str(totalD)
            if totalH > 0:
                SumName += ";O"
            if totalH > 1:
                SumName += str(totalH)
            SumName += ")"

            name = re.sub(
                "Alkyl_", "", self.lipidclass
            )  # to removed Alkyl labels from plasmalogens
            if (Nchains > 2) &  (self.lipidclass not in ["OAcylCeramide", "N_Acyl_PE", "N_Acyl_PS" ]):
                FullName = name + "(" + "_".join(ChainStrings) + ")"
            else:
                FullName = name + "(" + "/".join(ChainStrings) + ")"
            FullNameWithAdduct = FullName + " " + self.adduct

        else:  # for isoprene lipids, where N = # of isoprene units
            ChainStrings = str(self.chains[0][0])
            FullName = self.lipidclass + ChainStrings
            FullNameWithAdduct = FullName + " " + self.adduct

        atoms = self.MF()
        MASS = MW_list(atoms)
        FRAGMENTS = self.theoreticalDigest()
        FRAGMENTS.sort(key=lambda x: x[0], reverse=True)

        # to add together fragments that have same mass (e.g. two acyl groups of same composition)
        # if False:
        # 	i = 0
        # 	while i < (len(FRAGMENTS)-1):
        # 		if (abs(FRAGMENTS[i][0] - FRAGMENTS[i+1][0]) < 0.0001):  #combine if difference in mz is less than 0.1 millidaltons
        # 			FRAGMENTS[i][1] = FRAGMENTS[i][1] + FRAGMENTS[i+1][1]
        # 			FRAGMENTS[i][2] = FRAGMENTS[i][2] + "/" + FRAGMENTS[i+1][2]
        # 			del FRAGMENTS[i+1]
        # 		else:
        # 			i = i + 1

        # for MS3 case to collapse peaks in MS2 and separately collapse peaks in MS3 if the have the same MS2 precursor
        if True:
            FRAGMENTS_3 = []
            i = 0
            # pull out MS3 fragments and place in FRAGMENTS_3
            while i < (len(FRAGMENTS)):
                if re.match(".*ms3.*", FRAGMENTS[i][2]):
                    FRAGMENTS_3.append(FRAGMENTS[i])
                    del FRAGMENTS[i]
                else:
                    i = i + 1

            # collapse MS2 fragments
            i = 0
            while i < (len(FRAGMENTS) - 1):
                if (
                    abs(FRAGMENTS[i][0] - FRAGMENTS[i + 1][0]) < 0.0001
                ):  # combine if difference in mz is less than 0.1 millidaltons
                    FRAGMENTS[i][1] = FRAGMENTS[i][1] + FRAGMENTS[i + 1][1]
                    FRAGMENTS[i][2] = FRAGMENTS[i][2] + "/" + FRAGMENTS[i + 1][2]
                    del FRAGMENTS[i + 1]
                else:
                    i = i + 1

            # collapse MS3 fragments

            i = 0
            while i < (len(FRAGMENTS_3)):
                # Make 4th element (index 3) in list the MS2 precursor mass
                FRAGMENTS_3[i] = FRAGMENTS_3[i] + [
                    float(
                        re.sub(
                            "\{|\}",
                            "",
                            re.search(r"\{(.*?)\}", FRAGMENTS_3[i][2]).group(0),
                        )
                    )
                ]
                i = i + 1
            # sort by MS3 precursor mass
            FRAGMENTS_3.sort(key=lambda x: x[3], reverse=True)
            unique_MS2_prec = list(set([el[3] for el in FRAGMENTS_3]))

            # collapse only if we have the same precursor mass
            for mz in unique_MS2_prec:
                i = 0
                while i < (len(FRAGMENTS_3) - 1):
                    if FRAGMENTS_3[i][3] != mz or FRAGMENTS_3[i + 1][3] != mz:
                        i = i + 1
                        continue
                    elif (
                        abs(FRAGMENTS_3[i][0] - FRAGMENTS_3[i + 1][0]) < 0.0001
                    ):  # combine if difference in mz is less than 0.1 millidaltons
                        FRAGMENTS_3[i][1] = FRAGMENTS_3[i][1] + FRAGMENTS_3[i + 1][1]
                        FRAGMENTS_3[i][2] = (
                            FRAGMENTS_3[i][2] + "/" + FRAGMENTS_3[i + 1][2]
                        )
                        del FRAGMENTS_3[i + 1]
                    else:
                        i = i + 1

            # remove extra element of lists that had MS2 prec mass for tidyness
            i = 0
            while i < (len(FRAGMENTS_3)):
                del FRAGMENTS_3[i][3]
                i = i + 1

            # recombine MS2 and MS3 fragments
            FRAGMENTS = FRAGMENTS + FRAGMENTS_3
            # resort according to fragment mass
            FRAGMENTS.sort(key=lambda x: x[0], reverse=True)

        RECORD = []
        RECORD.append("Name: " + FullName + "\n")
        RECORD.append("Id: " + FullNameWithAdduct + "\n")
        RECORD.append("SumComposition: " + SumName + "\n")
        RECORD.append("MW: " + str(MASS) + "\n")
        RECORD.append("ExactMass: " + str(MASS) + "\n")
        RECORD.append(
            "PrecursorMz: "
            + str((MASS + ADDUCT[self.adduct]) / abs(ADDUCT_CHARGE[self.adduct]))
            + "\n"
        )
        RECORD.append("ADDUCT: " + self.adduct + "\n")
        RECORD.append("NCE: " + self.NCE + "\n")
        RECORD.append("CATEGORY: " + self.lipidclass + "\n")
        RECORD.append("CLASS: " + self.lipidclass + "\n")
        RECORD.append("LIPID MAPS CLASS: " + LIPID_MAPS[self.lipidclass] + "\n")
        RECORD.append("FORMULA: " + self.prettyFormula() + "\n")
        RECORD.append("NumPeaks: " + str(len(FRAGMENTS)) + "\n")

        for f in FRAGMENTS:
            RECORD.append(str(f[0]) + " " + str(f[1]) + " " + str(f[2]) + "\n")
        RECORD.append("\n\n")
        RECORDSTR = "".join(RECORD)

        return RECORDSTR

    def info(self):
        print(self.lipidclass, self.chains)
        atoms = self.MF()
        print("ATOMS:", atoms)
        print("MW:", MW_list(atoms))
        print("FORMULA:", self.prettyFormula())

    def theoreticalDigest(self):
        MASS = MW_list(self.MF())
        PREC = MASS + ADDUCT[self.adduct]

        FRAGMENTS = []
        # FRAGMENTS.append( [PREC, 1000, "pre"] )
        return FRAGMENTS

    def generateLibrary(self, target=None, mode="pos"):
        if target:
            handle = open(target, "a+")
        if mode == "pos":
            adduct_set = self.pos_adduct_set
        elif mode == "neg":
            adduct_set = self.neg_adduct_set
        # parent = self.__bases__[0]
        class_name = self.__class__.__name__
        for c in self.chain_sets:
            for adduct in adduct_set:
                self.set_chains_and_adduct(class_name, c, adduct=adduct)
                content = self.printNist()
                if target:
                    handle.write(content)
                else:
                    sys.stdout.write(content)
        if target:
            handle.close()


# DG, PA, PC, PE, PG, PI, PS, BMP, HexDG, MMPE, DMPE, DG, FAFHFA
class GPL(Lipid):
    chain1_ranges = []
    for c in [2] + list(range(10, 22)) + list(range(22, 34, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain1_ranges.append([c, d, h])

    chain2_ranges = []
    for c in [2] + list(range(10, 22)) + list(range(22, 34, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain2_ranges.append([c, d, h])

    chain_sets = []
    for c1 in chain1_ranges:
        for c2 in chain2_ranges:
            # if c2[0] < c1[0]: continue
            # if (c2[0] == c1[0] and c2[1] < c1[0]): continue
            chain_sets.append([c1, c2])


# Alkyl_PC, Alkyl_PE, Alkyl_PS
class alkylGPL(Lipid):
    chain1_ranges = []
    for c in [2] + list(range(10, 23)):
        for d in range(0, 3):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain1_ranges.append([c, d, h])

    chain2_ranges = []
    for c in [2] + list(range(10, 22)) + list(range(22, 34, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain2_ranges.append([c, d, h])

    chain_sets = []
    for c1 in chain1_ranges:
        for c2 in chain2_ranges:
            chain_sets.append([c1, c2])


# TG,
class Triglyceride(Lipid):
    chain1_ranges = []
    for c in [2] + list(range(10, 26, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain1_ranges.append([c, d, h])

    chain2_ranges = chain1_ranges
    chain3_ranges = chain1_ranges

    chain_sets = []
    for c1 in chain1_ranges:
        for c2 in chain2_ranges:
            for c3 in chain3_ranges:
                if c2[0] < c1[0]: continue
                if (c2[0] == c1[0] and c2[1] < c1[1]): continue
                if c3[0] < c2[0]: continue
                if (c3[0] == c2[0] and c3[1] < c2[1]): continue
                chain_sets.append([c1, c2, c3])


# N_Acyl_PE, N_Acyl_PS
class NAcylGPL(Lipid):
    chain1_ranges = []
    for c in [2] + list(range(10, 26, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain1_ranges.append([c, d, h])

    chain2_ranges = []
    for c in [2] + list(range(10, 26, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain2_ranges.append([c, d, h])

    chain3_ranges = []
    for c in [2] + list(range(10, 26, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain3_ranges.append([c, d, h])

    chain_sets = []
    for c1 in chain1_ranges:
        for c2 in chain2_ranges:
            for c3 in chain3_ranges:
                if c2[0] < c1[0]:
                    continue
                if c2[0] == c1[0] and c2[1] < c1[0]:
                    continue
                chain_sets.append([c1, c2, c3])


# AcGM2, AcGM3, GcGM2, GcGM3, MIP2C, MIPC, HexCer, LacCer, PI_Cer, Sulfatide, SM, Ceramide, Ceramide_P
# Note that for Sphingolipids, the hydroxyl "head group" does not count as a hydroxyl for the purposes of this script


class SphingoLipid(Lipid):
    # chain1_ranges = []
    # for c in (range(14,23)):
    # 	for d in range(0,3):
    # 		for h in range(0,3):
    # 			chain1_ranges.append([c,d,h])

    #
    # chain2_ranges = []
    # for c in ([2] + range(10,22) + range(22,34,2)):
    # 	for d in range(0,7):
    # 		if (c > 5 and c < 22 and d > (c-5)/3) or (c < 6 and d > 0): continue
    # 		for h in range(0,2):
    # 			chain2_ranges.append([c,d,h])
    #
    #
    chain1_ranges = []
    for c in range(14, 23):
        for d in range(0, 2):
            for h in range(0, 3):
                chain1_ranges.append([c, d, h])

    chain2_ranges = []
    for c in [2] + list(range(14, 25)):
        for d in range(0, 7):
            for h in range(0, 2):
                chain2_ranges.append([c, d, h])

    chain_sets = []
    for c1 in chain1_ranges:
        for c2 in chain2_ranges:
            chain_sets.append([c1, c2])


class AcylSphingoLipid(Lipid):
    chain1_ranges = []
    for c in range(14, 21):
        for d in range(0, 2):
            for h in range(1, 2):
                chain1_ranges.append([c, d, h])

    chain2_ranges = []
    for c in [2] + list(range(14, 25)):
        for d in range(0, 7):
            for h in range(0, 2):
                chain2_ranges.append([c, d, h])

    chain3_ranges = []
    for c in [2] + list(range(14, 25)):
        for d in range(0, 7):
            for h in range(0, 2):
                chain3_ranges.append([c, d, h])

    chain_sets = []
    for c1 in chain1_ranges:
        for c2 in chain2_ranges:
            for c3 in chain3_ranges:
                chain_sets.append([c1, c2, c3])

    #chain_sets = [[[18,1,1], [17,0,0],[18,1,0]]]


# LPA, LPC, LPE, LPG, LPI, LPS
class LysoGPL(Lipid):
    chain1_ranges = []
    for c in list(range(2, 22)) + list(range(22, 34, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain1_ranges.append([c, d, h])
    chain_sets = []
    for c1 in chain1_ranges:
        chain_sets.append([c1, [0, 0, 0]])
        chain_sets.append([[0, 0, 0], c1])


# MG
class MAG(Lipid):
    chain1_ranges = []
    for c in list(range(2, 22)) + list(range(22, 34, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain1_ranges.append([c, d, h])
    chain_sets = []
    for c1 in chain1_ranges:
        chain_sets.append([c1])


# Alkyl_LPC, ALkyl_LPE, Alkyl_LPS
class AlkylLysoGPL(Lipid):
    chain1_ranges = []
    for c in [2] + list(range(10, 23)):
        for d in range(0, 3):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain1_ranges.append([c, d, h])
    chain_sets = []
    for c1 in chain1_ranges:
        chain_sets.append([c1])


# LysoHexCer, LysoSM, LCB
# Note that for Sphingolipids, the hydroxyl "head group" does not count as a hydroxyl for the purposes of this script
class LysoSphingoLipid(Lipid):
    chain1_ranges = []
    for c in range(14, 23):
        for d in range(0, 3):
            for h in range(0, 3):
                chain1_ranges.append([c, d, h])
    chain_sets = []
    for c1 in chain1_ranges:
        chain_sets.append([c1])


# CE, ErgE, Ethanolamine, FA, Carn,
class singleAcyl(Lipid):
    chain1_ranges = []
    for c in list(range(2, 22)) + list(range(22, 34, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0, 1]:
                chain1_ranges.append([c, d, h])
    chain_sets = []
    for c1 in chain1_ranges:
        chain_sets.append([c1])


# CL, LysoCL
class CardioLipin(Lipid):
    chain1_ranges = []
    for c in [2] + list(range(14, 26, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain1_ranges.append([c, d, h])
    chain2_ranges = chain1_ranges
    chain3_ranges = chain1_ranges
    chain4_ranges = chain1_ranges

    chain_sets = []
    for c1 in chain1_ranges:
        for c2 in chain2_ranges:
            for c3 in chain3_ranges:
                for c4 in chain4_ranges:
                    if c2[0] < c1[0]: continue
                    if (c2[0] == c1[0] and c2[1] < c1[1]): continue
                    if c3[0] < c2[0]: continue
                    if (c3[0] == c2[0] and c3[1] < c2[1]): continue
                    if c3[0] < c2[0]: continue
                    if (c3[0] == c2[0] and c3[1] < c2[1]): continue
                    if c4[0] < c3[0]: continue
                    if (c4[0] == c3[0] and c3[1] < c2[1]): continue

                    # if c2[0] < c1[0]: continue
                    # if (c2[0] == c1[0] and c2[1] < c1[0]): continue
                    # if c4[0] < c3[0]: continue
                    # if (c4[0] == c3[0] and c3[1] < c1[0]): continue
                    # if (c3[0] < c1[0] or (c3[0] == c1[0] and c3[1] < c3[1])): continue # use for Maven, but not database
                    chain_sets.append([c1, c2, c3, c4])


class LysoCardioLipin(Lipid):
    chain1_ranges = []
    for c in [2] + list(range(14, 26, 2)):
        for d in range(0, 7):
            if (c > 5 and c < 22 and d > (c - 5) / 3) or (c < 6 and d > 0):
                continue
            for h in [0]:
                chain1_ranges.append([c, d, h])
    chain2_ranges = chain1_ranges
    chain3_ranges = chain1_ranges

    chain_sets = []
    for c1 in chain1_ranges:
        for c2 in chain2_ranges:
            for c3 in chain3_ranges:
                if c2[0] < c1[0]: continue
                if (c2[0] == c1[0] and c2[1] < c1[1]): continue
                if c3[0] < c2[0]: continue
                if (c3[0] == c2[0] and c3[1] < c2[1]): continue
                if c3[0] < c2[0]: continue
                # if c2[0] < c1[0]: continue
                # if (c2[0] == c1[0] and c2[1] < c1[0]): continue
                # if c3[0] < c2[0]: continue
                # if (c3[0] == c2[0] and c3[1] < c1[0]): continue
                chain_sets.append([c1, c2, c3])


class Isoprenoid(Lipid):
    chain1_ranges = []
    for c in range(1, 21):
        chain1_ranges.append([c, 0, 0])
    chain_sets = []
    for c1 in chain1_ranges:
        chain_sets.append([c1])
