import datetime
import argparse
import os
import pandas as pd
from nltk import rte_features

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

    parser.add_argument(
        "-r",
        "--retention-times",
        help="retention time file name",
        dest="retention_time_file_name",
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
        print("-r <retention times>:")
        print("\tSupply path and filename for a .tsv containing known retention times\n")
        #Needs Name RT or Name and RT+rt_min+rt_max or Name and rt_min_rt_max?
        print("\tFile must contain a column with \"name\" and \"rt\" and/or\n")
        print("\t\"name\" and both of \"rt_max\" and \"rt_min\"\n")
        print("\tDEFAULT: No retention times will be associated.\n")

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



    if args.retention_time_file_name != "":
        if not os.path.exists(args.retention_time_file_name):
            print("Supplied retention time files does not exist. Exiting program.")
            exit(1)

        retention_times_df = pd.read_csv(args.retention_time_file_name, delimiter = '\t')
        retention_times_df.columns = retention_times_df.columns.str.lower() #make column names lower case
        rt_col_names = retention_times_df.columns.tolist()

        #check that we have Name and RT or both of RT_min and RT_max (or all three)
        if ("name" not in rt_col_names) and (("rt" in rt_col_names) or ("rt_max" not in rt_col_names and "rt_min" not in rt_col_names)):
            print("Supplied retention time file is improperly formatted.\n")
            print("\tSupply path and filename for a .tsv containing known retention times\n")
            print("\tFile must contain a column with \"name\" and \"rt\" and/or\n")
            print("\tboth of \"rt_max\" and \"rt_min\"\n")
            print("\tExiting program.")
            exit(1)
        else:
            retention_times_df = retention_times_df.loc[:, retention_times_df.columns.isin(["name", "rt", "rt_min", "rt_max"])]
    else:
        retention_times_df = None

    lipid_classes = []
    if len(args.lipid_classes) == 0:
        lipid_classes = ALL_LIPID_CLASSES.copy()
    else:
        lipid_classes_raw = [
            re.sub(r'["\']', "", x.strip()) for x in args.lipid_classes.split(",")
        ]
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

    lipid_classes_list = list(lipid_classes)
    lipid_classes_list.sort() #put lipid classes in alphabetical order
    for lipid_class in lipid_classes_list:
        msg = "Enumerating spectra for " + lipid_class + "  ... "
        instantiate_str = lipid_class + " = " + lipid_class + "()"
        sys.stdout.write(msg)
        eval_str = (f"{lipid_class}.generateLibrary(target='{output_msp_file}', mode='{ion}', retention_times=retention_times_df)"
                    )


        exec(instantiate_str)
        eval(eval_str)
        sys.stdout.write(" DONE\n")

    print("All Processes Successfully Completed!")
    exit(0)
