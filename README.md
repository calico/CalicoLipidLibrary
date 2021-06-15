# CalicoLipidLibrary
Scripts to create in silico lipid fragmentation spectra, based on analysis of chemical standards and relevant scientific literature.

# Installation

1.  clone this repository.

2. `cd CalicoLipidLibrary`

3. `python setup.py sdist bdist_wheel`

# Generate A Lipid Library MSP File

This is the primary usage of this package.  Produce an msp file for all desired in the appropriate ionization mode:
```
# all negative mode lipids
python2.7 generateDB.py -m "neg"

# PC in positive mode
python2.7 generateDB.py -m "pos" -c "PC"

# PC, PE, PS in positive mode. ImaginaryClass is silently removed
python2.7 generateDB.py -m "pos" -c "PC,PE,ImaginaryClass,PS"

# Sphingolipids neg mode, into file named 'neg_sphingolipids.msp'
python2.7 generateDB.py -m "neg" -c "AcGM2, AcGM3, GcGM2, GcGM3, MIP2C, MIPC, HexCer, LacCer, Sulfatide, SM, Ceramide, Ceramide_P" -n "neg_sphingolipids"

#  Glycophospholipids pos mode, into file named "glycophospholipids_pos.msp", in directory "/Users/SomeUser/calicolipidlibrary_output"
python2.7 generateDB.py -m "pos" -o "/Users/SomeUser/calicolipidlibrary_output" -n "glycophospholipids_pos.msp" -c "DG, PA, PC, PE, PG, PI, PS, BMP, HexDG, MMPE, DMPE, DG, FAFHFA"
```

You may also run the script with no arguments to print a usage message.
```
# Print usage information
python2.7 generateDB.py
```
```
Usage:

python generateDB.py -m <pos|neg> -o <output-path> -c <lipid-classes>

-m <pos|neg>:
	enumerate positive or negative mode spectra.
	DEFAULT: neg (negative ionization mode)
-o <output-path>:
	Supply desired output path.
	DEFAULT: Current directory
-n <output-file-name>:
	Name of output lipids library file.
	DEFAULT: <date>-Calico-Lipids-<all|<classes>>-<ion>.msp
-c <lipid-classes>:
	Supply desired lipid classes, in a comma-delimited string (e.g., 'PS,PC,PE'
	DEFAULT: all available lipid classes.

AVAILABLE LIPID CLASSES:
CPI
DGTS
DG
MMPE
Taurine
CPE
Sulfatide
GcGM3
FA
LysoCL
MGDG
AcGM3
AcGM2
PS
LacCer
N_Acyl_PE
PC
MIP2C
LysoHexCer
ErgE
BMP
N_Acyl_PS
Carn
PG
MIPC
DMPE
TG
PI
FAHFA
CDP_DG
Ceramide_P
LysoSM
CL
LCB_P
PA
PE
LPS
CE
BDP
LPI
LPG
LPE
Ethanolamine
LPC
Ceramide
LPA
MG
Alkyl_PC
LysoCPI
Alkyl_PE
GcGM2
LysoCPE
DGDG
Alkyl_LPS
HemiBMP
Alkyl_PS
LCB
SM
Alkyl_LPE
HexCer
Alkyl_LPC
```
# Print Single Spectrum

Running python in an interactive shell, you may also print the theoretical MS/MS spectrum of a single lipid class+adduct:

1. Launch interactive shell. Currently, only python 2 is supported.
```
python2.7
```

2. Import the library into the python session.
```
import calicolipidlibrary
```

3. Print a spectrum of a single lipid class/adduct form of interest.

If the lipid class or adduct form is not supported, a list of all supported classes and adducts will appear.
If the fragmentation does not make sense for the given chain combination and structural composition of the lipid, an error will display.

Usage:
```
calicolipidlibrary.print_spectrum(<class_name>, <acyl_chains>, <adduct_name>)

<class_name>: Supported lipid class. See examples below for complete list.

<acyl_chains>: [<chain_1>, <chain_2>, <chain_3>, ...]
within <chain_1>: [<num_carbons>, <num_double_bonds>, <num oxidations>]
e.g. [[18,1,1], [16:0,2], [22,6,0]]
sn1 is 18:1 with one oxidation (18:1;O)
sn2 is 16:0 with two oxidations (16:0;O2)
sn3 is 22:6 with no oxidations (22:6)

Note that the p- forms of Alkyl_PC and Alkyl_PE imply an additional double bond,
and that when supplying oxidations to sphingolipids, an additional OH will be added from the hydroxyl head group.

<adduct_name>: Supported adduct form. valid adduct forms:
"[M+H]+", "[M+Na]+",  "[M+NH4]+", "[M+K]+", "[M+2H]2+", "[M+Li]+"
"[M-H]-", "[M+Cl]-","[M+37Cl]-", "[M+FA-H]-", "[M-CH3]-", "[M+AcOH-H]-", "[M-2H]2-", "[M-3H]3-"

Not all adduct names are valid for all lipid classes.
```
Here is a list of valid `print_spectrum()` commands, showing one valid spectrum per class:
```
calicolipidlibrary.print_spectrum("AcGM2", [[18,1,1], [18,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("AcGM3", [[18,1,1], [24,1,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("Alkyl_LPC", [[18,1,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("Alkyl_LPE", [[18,1,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("Alkyl_LPS",[[16,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("Alkyl_PC", [[16,1,0], [18,0,0]], "[M+FA-H]-")
calicolipidlibrary.print_spectrum("Alkyl_PE", [[18,1,0], [22,6,0]], "[M+AcOH-H]-")
calicolipidlibrary.print_spectrum("Alkyl_PS", [[16,0,0], [18,2,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("BDP", [[18,1,0], [18,0,0], [16,1,0],[16,0,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("BMP", [[18,1,0], [18,3,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("Carn", [[16,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("CDP_DG", [[18,1,0], [18,2,0]], "[M+Na]+")
calicolipidlibrary.print_spectrum("CE", [[18,2,1]],"[M+Na]+")
calicolipidlibrary.print_spectrum("Ceramide",[[18,1,1], [16,0,1]], "[M+H]+")
calicolipidlibrary.print_spectrum("Ceramide_P", [[18,1,1], [12,0,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("CL", [[16,0,0], [18,1,0], [16,1,0], [18,0,0]], "[M-2H]2-")
calicolipidlibrary.print_spectrum("CPE", [[18,1,1], [24,1,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("CPI", [[18,1,1], [24,0,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("DG", [[16,0,0], [20,1,0]], "[M+NH4]+")
calicolipidlibrary.print_spectrum("DGDG", [[18,3,0], [16,3,0]], "[M+NH4]+")
calicolipidlibrary.print_spectrum("DGTS", [[16,0,0], [16,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("DMPE", [[16,0,0], [16,1,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("ErgE", [[16,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("Ethanolamine", [[16,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("FA", [[16,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("FAHFA", [[18,0,0], [18,1,0]], "[M+Na]+")
calicolipidlibrary.print_spectrum("GcGM2", [[18,1,1], [24,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("GcGM3", [[18,1,1], [24,1,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("HemiBMP", [[18,0,0], [22,6,0], [2,0,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("HexCer", [[18,1,2], [18,2,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("LacCer", [[18,1,1], [24,1,0]], "[M+FA-H]-")
calicolipidlibrary.print_spectrum("LCB", [[18,1,1]], "[M+H]+")
calicolipidlibrary.print_spectrum("LCB_P", [[20,1,1]], "[M+H]+")
calicolipidlibrary.print_spectrum("LPA", [[16,0,0],[0,0,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("LPC", [[18,0,0],[0,0,0]], "[M+Na]+")
calicolipidlibrary.print_spectrum("LPE", [[18,0,0],[0,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("LPG", [[18,1,0],[0,0,0]], "[M+Na]+")
calicolipidlibrary.print_spectrum("LPI", [[16,0,0],[0,0,0]], "[M+Na]+")
calicolipidlibrary.print_spectrum("LPS", [[18,1,0],[0,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("LysoCL", [[16,0,0], [18,1,0], [18,2,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("LysoCPE", [[18,1,1]], "[M-H]-")
calicolipidlibrary.print_spectrum("LysoCPI", [[18,1,1]], "[M-H]-")
calicolipidlibrary.print_spectrum("LysoHexCer", [[16,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("LysoSM",[[16,0,0]], "[M+Na]+")
calicolipidlibrary.print_spectrum("MG", [[14,0,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("MGDG", [[18,3,0], [16,3,0]], "[M+NH4]+")
calicolipidlibrary.print_spectrum("MIP2C", [[18,1,1], [24,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("MIPC", [[18,0,2], [26,0,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("MMPE", [[16,1,0], [16,1,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("N_Acyl_PE", [[16,1,0], [18,2,0], [20,2,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("N_Acyl_PS", [[16,0,0], [18,1,0], [18,2,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("PA",[[16,0,0], [18,1,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("PC",[[18,0,0], [20,4,0]], "[M+FA-H]-")
calicolipidlibrary.print_spectrum("PE",[[16,0,0], [18,1,0]], "[M-H]-")
calicolipidlibrary.print_spectrum("PG",[[16,0,0], [18,1,0]], "[M+Na]+")
calicolipidlibrary.print_spectrum("PI", [[20,4,0], [16,0,0]], "[M+Na]+")
calicolipidlibrary.print_spectrum("PS",[[16,0,0], [18,1,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("SM",[[18,1,0], [20,4,3]], "[M+Na]+")
calicolipidlibrary.print_spectrum("Sulfatide", [[18,0,1], [16,1,0]], "[M+FA-H]-")
calicolipidlibrary.print_spectrum("Taurine", [[18,0,0]], "[M+H]+")
calicolipidlibrary.print_spectrum("TG", [[16,0,0], [18,1,0], [18,0,0]], "[M+NH4]+")
```