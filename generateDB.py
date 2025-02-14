import datetime
import argparse
import os
from calicolipidlibrary import *


def cli_parser():
    parser = argparse.ArgumentParser(
        description="entry point to create calico lipid libraries",
        usage="python generateDB.py -m [pos|neg] -o [output_directory",
    )

    parser.add_argument(
        "-m", "--ion-mode", help="ion_mode", dest="ion_mode", default="neg"
    )

    parser.add_argument(
        "-o",
        "--output-path",
        help="output folder",
        dest="output_folder",
        default=os.getcwd(),
    )

    parser.add_argument(
        "-c",
        "--lipid-classes",
        help="list of lipid classes to include",
        dest="lipid_classes",
        default="",
    )

    parser.add_argument(
        "-n",
        "--output-name",
        help="output file name",
        dest="output_file_name",
        default="",
    )

    return parser


if __name__ == "__main__":
    try:
        arg1 = sys.argv[1]
    except IndexError:
        print("Usage:\n")
        print("python generateDB.py -m <pos|neg> -o <output-path> -c <lipid-classes>\n")
        print("-m <pos|neg>:")
        print("\tenumerate positive or negative mode spectra.")
        print("\tDEFAULT: neg (negative ionization mode)")
        print("-o <output-path>:")
        print("\tSupply desired output path.")
        print("\tDEFAULT: Use this directory")
        print("-n <output-file-name>:")
        print("\tName of output lipids library file.")
        print("\tDEFAULT: <date>-Calico-Lipids-[all|<lipid-classes>]-<ion>.msp")
        print("-c <lipid-classes>:")
        print(
            "\tSupply desired lipid classes, in a comma-delimited string (e.g., 'PS,PC,PE'"
        )
        print("\tDEFAULT: all available lipid classes.\n")
        print("AVAILABLE LIPID CLASSES:")
        for lipid_class in ALL_LIPID_CLASSES:
            print(lipid_class)
        exit(0)

    args = cli_parser().parse_args()

    if not os.path.exists(args.output_folder):
        print("Supplied output folder path does not exist. Exiting program.")
        exit(1)

    path = args.output_folder

    ion = None
    if args.ion_mode == "neg" or args.ion_mode == "pos":
        ion = args.ion_mode
    else:
        print(
            "Supplied ion mode argument could not be interpreted,\n"
            "please use 'neg' for negative mode or 'pos' for positive mode.\n"
            "Exiting program."
        )

        exit(1)

    lipid_classes = []
    if len(args.lipid_classes) == 0:
        lipid_classes = ALL_LIPID_CLASSES.copy()
    else:
        lipid_classes_raw = [x.strip() for x in args.lipid_classes.split(",")]
        for lipid_class in lipid_classes_raw:
            if lipid_class in ALL_LIPID_CLASSES:
                lipid_classes.append(lipid_class)

    file_name = args.output_file_name  # file name only, no system path
    if len(file_name) == 0:
        now = datetime.datetime.now()
        date = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)
        fnbase = "-Calico-Lipids-"
        class_str = "all"

        if len(args.lipid_classes) != 0:
            class_str = "-".join(lipid_classes)

        file_name = date + fnbase + class_str + "-" + ion + ".msp"

    if not file_name.endswith(".msp"):
        file_name = file_name + ".msp"

    output_msp_file = os.path.join(path, file_name)  # full system path

    print("====================================")
    print("PARAMETERS:")
    print("output file: " + output_msp_file)
    print("ion mode: " + ion)
    if len(args.lipid_classes) == 0:
        print("lipid classes: ALL")
    else:
        print("lipid classes:")
        for lipid_class in lipid_classes:
            print(lipid_class)
    print("====================================")

    open(output_msp_file, "w")

    for lipid_class in lipid_classes:
        msg = "Enumerating spectra for " + lipid_class + "  ... "
        sys.stdout.write(msg)
        eval_str = (
            lipid_class
            + ".generateLibrary(target='"
            + output_msp_file
            + "', mode='"
            + ion
            + "')"
        )
        eval(eval_str)
        sys.stdout.write(" DONE\n")

    print("All Processes Successfully Completed!")
    exit(0)
