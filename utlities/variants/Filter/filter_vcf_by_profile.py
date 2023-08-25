#!/usr/bin/env python3
import os
import sys
import yaml
import gzip
from icecream import ic

if len(sys.argv) < 3:
    print(f"this script is to filter VCF by a configure file")
    sys.exit(f"python3 {sys.argv[0]} vcf profile")

def load_config(profile_file):
    # Read YAML file
    with open(profile_file, 'r') as stream:
        try:
            data_loaded = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            sys.exit(exc)
    return data_loaded

def genotyping(content):
    """
    0/0. 0/1, 1/1, 1/2, 2/2, 1|1, 2|2 ...
    """
    return content[0:3]

def smart_open(file, opera):
    if opera == 'r':
        if os.path.exists(file) ==False:
            print("Can not open file {}".format(file))
            exit()
        else:
            if file.endswith(".gz"):
                out = gzip.open(file, 'rt')
            else:
                out = open(file, 'r')
    elif opera == 'w':
        if file.endswith(".gz"):
            out = gzip.open(file, 'wt')
        else:
            out = open(file, 'w')
    return out


def main():
    sample_order = {} # link the sample order to columns
    config = load_config(sys.argv[2])
    filtering = config['filter'] # yes means True, no means False in YAML
    ic(filtering)
    case = config['case']
    case_number = len(case)
    positive_in_case = config['positive_in_case']
    control = config['control']
    control_number = len(control)
    missing_rate = config['missing']

    with smart_open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            if i.startswith("##"):
                continue
            elif filtering and tmp[6] == 'Filter':
                continue
            elif i[0] == "#": # this is header, get sample order from this line
                for order in range(9,len(tmp)):
                    #ic(order, tmp[order])
                    sample_order[order] = tmp[order]
            else:
                # content
                case_11 = 0
                control_00 = 0
                control_11 = 0
                control_22 = 0
                missing = 0
                for order,content in enumerate(tmp[9:], start=9):
                    genotype = genotyping(content)
                    sample = sample_order[order]
                    #ic(order, sample)
                    if sample in case:
                        if genotype == '1/1' or genotype == '1|1':
                            case_11 += 1
                    else:
                        if genotype == '0/0' or genotype == '0|0':
                            control_00 += 1
                        elif genotype == '1/1' or genotype == '1|1':
                            control_11 += 1
                        elif genotype == '2/2' or genotype == '2|2':
                            control_22 += 1
                        elif genotype == './.':
                            missing += 1
                # alt in case; all samples in control are 0/0 or missing
                if case_11 >= positive_in_case and control_22 == 0 and missing <= control_number * missing_rate and control_00 + missing == control_number:
                    #print("case: 1/1; control 0/0")
                    print(i.strip())
                # alt in case; all samples in control are 2/2 or missing
                elif case_11 >= positive_in_case and control_11 == 0 and missing <= control_number * missing_rate and control_22 + missing == control_number:
                    #print("case: 1/1; control 2/2")
                    print(i.strip())

if __name__ == '__main__':
    main()
