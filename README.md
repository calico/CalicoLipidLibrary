# CalicoLipidLibrary
Scripts to create in silico lipid fragmentation spectra, based on analysis of chemical standards and relevant scientific literature.

# Installation

1.  clone this repository.

2. `cd CalicoLipidLibrary`

3. `python setup.py sdist bdist_wheel`

# Generate Entire Lipid Library

This is the primary usage of this package.  Produce an msp file for all desired/classes and adducts, like so:

`TODO`

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

TODO

calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("",,"")
calicolipidlibrary.print_spectrum("Ceramide",[[18,1,1], [16,0,0]], "[M+H]+")
```