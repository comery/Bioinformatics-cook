#!/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/bin/python
import os
import sys
import time
import gzip
import argparse
import subprocess
from bold_identification.BOLD_identification import (
    main as bold_identification,
)

t = time.time()

try:
    import Bio
except:
    sys.exit("package biopython not found! Please install it!")
else:
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna


###############################################################################
#####------------------------- parameters --------------------------------#####

## common group  ##
common_parser = argparse.ArgumentParser(add_help=False)

common_group = common_parser.add_argument_group("common arguments")

common_group.add_argument(
    "-outpre",
    metavar="<STR>",
    required=True,
    help="prefix for output files",
)
## index group ##
index_parser = argparse.ArgumentParser(add_help=False)

index_group = index_parser.add_argument_group("index arguments")

index_group.add_argument(
    "-index",
    metavar="INT",
    type=int,
    required=True,
    help="the length of tag sequence in the ends of primers",
)

# Software group
soft_parser = argparse.ArgumentParser(add_help=False)

soft_group = soft_parser.add_argument_group("software path")

soft_group.add_argument(
    "-vsearch",
    metavar="<STR>",
    help="vsearch path" + "(only needed if vsearch is not in $PATH)",
)

soft_group.add_argument(
    "-threads", metavar="<INT>", default=2, help="threads for vsearch, default=2"
)

soft_group.add_argument(
    "-cid",
    metavar="FLOAT",
    type=float,
    default=0.98,
    dest="cluster_identity",
    help="identity for clustering, default=0.98",
)

## filter group  ##
filter_parser = argparse.ArgumentParser(
    add_help=False,
    description="Use the whole raw" + "dataset (Only adapters should be removed)!",
)

filter_group = filter_parser.add_argument_group("filter arguments")

filter_group.add_argument(
    "-raw",
    metavar="<STR>",
    required=True,
    help="input raw Single-End fastq file, and"
    + " only\nadapters should be removed; supposed on\n"
    + "Phred33 score system (BGISEQ-500)",
)

filter_group.add_argument(
    "-phred",
    metavar="<INT>",
    type=int,
    dest="phred",
    choices=[33, 64],
    default=33,
    help="Phred score system, 33 or 64, default=33"
)

filter_group.add_argument(
    "-e",
    metavar="<INT>",
    type=int,
    dest="expected_err",
    help="expected error threshod, default=10\n"
    + "see more: http://drive5.com/usearch/manual/exp_errs.html",
)

filter_group.add_argument(
    "-q",
    metavar="<INT>",
    type=int,
    dest="quality",
    nargs=2,
    help="filter by base quality; for example: '20 5' means\n"
    + "dropping read which contains more than 5 percent of \n"
    + "quality score < 20 bases.",
)

filter_group.add_argument(
    "-trim",
    dest="trim",
    action="store_true",
    help="whether to trim 5' end of read, it adapts to -e mode\n"
    + "or -q mode",
)
filter_group.add_argument(
    "-n",
    metavar="<INT>",
    type=int,
    default=1,
    help="remove reads containing [INT] Ns, default=1",
)

# ------------------------------------------------------------------------------

## assign group ##
assign_parser = argparse.ArgumentParser(
    add_help=False,
    description="assing clean reads to"
    + "samples by unique tag sequence"
    +"with 100% similarity",
)

assign_group = assign_parser.add_argument_group("assign arguments")

assign_group.add_argument(
    "-primer",
    metavar="<STR>",
    required=True,
    help="taged-primer list, on following format:\n"
    + "Rev001   AAGCTAAACTTCAGGGTGACCAAAAAATCA\n"
    + "For001   AAGCGGTCAACAAATCATAAAGATATTGG\n"
    + "...\n"
    + "this format is necessary!",
)

assign_group.add_argument(
    "-outdir",
    metavar="<STR>",
    default="assigned",
    help="output directory for assignment," + 'default="assigned"',
)

assign_group.add_argument(
    "-tmis",
    metavar="<INT>",
    type=int,
    dest="tag_mismatch",
    default=0,
    help="mismatch number in tag when demultiplexing, default=0",
)

assign_group.add_argument(
    "-pmis",
    metavar="<INT>",
    type=int,
    dest="primer_mismatch",
    default=1,
    help="mismatch number in primer when demultiplexing, default=1",
)

## only assign need
only_assign_parser = argparse.ArgumentParser(add_help=False)
only_assign_group = only_assign_parser.add_argument_group(
    "when only run assign arguments"
)

only_assign_group.add_argument(
    "-fq", metavar="<STR>", required=True, help="cleaned fastq file (*.fq.gz, *.fq)"
)

# ------------------------------------------------------------------------------

## assembly group ##
assembly_parser = argparse.ArgumentParser(
    description="Due to connect HIFI barcode sequence by overlaping two"
    + "consensus sequences which generated from clustering method,"
    + "in this program I use VSEARCH.  You can define the length of overlap,"
    + "how many reads used to make clusters, and whether to check codon"
    + "translation for PCG.",
    add_help=False,
)

assembly_group = assembly_parser.add_argument_group("assembly arguments")

assembly_group.add_argument(
    "-min",
    metavar="INT",
    type=int,
    default=80,
    dest="min_overlap",
    help="minimun length of overlap, default=80",
)

assembly_group.add_argument(
    "-max",
    metavar="INT",
    type=int,
    default=90,
    dest="max_overlap",
    help="maximum length of overlap, default=90",
)

assembly_group.add_argument(
    "-oid",
    metavar="FLOAT",
    type=float,
    default=0.95,
    dest="overlap_identity",
    help="minimun similarity of overlap region, default=0.95",
)

assembly_group.add_argument(
    "-tp",
    metavar="INT",
    type=int,
    dest="cluster_number_needKeep",
    help="how many clusters will be used in" + "assembly, recommend 2",
)

assembly_group.add_argument(
    "-ab",
    metavar="INT",
    type=int,
    dest="abundance_threshod",
    help="keep clusters to assembly if its abundance >=INT ",
)

assembly_group.add_argument(
    "-seqs_lim",
    metavar="INT",
    type=int,
    default=0,
    help="reads number limitation. by default,\n" + "no limitation for input reads",
)

assembly_group.add_argument(
    "-len",
    metavar="INT",
    type=int,
    default=400,
    dest="standard_length",
    help="standard read length, default=400",
)

assembly_group.add_argument(
    "-ds",
    dest="drop_short_read",
    action="store_true",
    help="drop short reads away before assembly",
)

assembly_group.add_argument(
    "-mode",
    metavar="INT",
    type=int,
    choices=[1, 2],
    default=1,
    help="1 or 2; modle 1 is to cluster and keep\n"
    + "most [-tp] abundance clusters, or clusters\n"
    + "abundance more than [-ab], and then make a \n"
    + "consensus sequence for each cluster.\n"
    + "modle 2 is directly to make only one consensus\n"
    + "sequence without clustering. default=1",
)

assembly_group.add_argument(
    "-rc",
    dest="reads_check",
    action="store_true",
    help="whether to check amino acid\n"
    +"translation for reads, default not",
)

# translation need
trans_parser = argparse.ArgumentParser(add_help=False)
trans_group = trans_parser.add_argument_group(
    "translation arguments(when set -rc or -cc)"
)

trans_group.add_argument(
    "-codon",
    metavar="INT",
    type=int,
    dest="codon_table",
    default=5,
    help="codon usage table used to check" + "translation, default=5",
)

trans_group.add_argument(
    "-frame",
    metavar="INT",
    type=int,
    choices=[0, 1, 2],
    default=1,
    help="start codon shift for amino acid" + "translation, default=1",
)

## only assembly need
only_assembly_parser = argparse.ArgumentParser(add_help=False)
only_assembly_group = only_assembly_parser.add_argument_group(
    "only run assembly arguments(not all)"
)

only_assembly_group.add_argument(
    "-list",
    metavar="FILE",
    type=str,
    required=True,
    help="input file, fastq file list. [required]",
)

polish_parser = argparse.ArgumentParser(
    description="polish all assemblies, \n"
    + "to make a confident COI barcode"
    + "reference.",
    add_help=False,
)
polish_group = polish_parser.add_argument_group(
    "polish arguments"
)

polish_group.add_argument(
    "-i",
    metavar="STR",
    type=str,
    dest="coi_input",
    required=True,
    help="COI barcode assemblies (fasta)",
)

polish_group.add_argument(
    "-o",
    metavar="STR",
    type=str,
    dest="coi_output",
    help="polished COI barcode assemblies (fasta)",
)

polish_group.add_argument(
    "-tag",
    metavar="STR",
    type=str,
    dest="sampleMark",
    help="add a mark for each sampel, like: >MARK_001;xxx",
)

polish_group.add_argument(
    "-cc",
    dest="coi_check",
    action="store_false",
    help="whether to check final COI contig's\n"
    + "amino acid translation, default yes",
)

polish_group.add_argument(
    "-cov",
    metavar="INT",
    type=int,
    dest="min_coverage",
    default=5,
    help="minimun coverage of 5' or 3' end allowed, default=5",
)

polish_group.add_argument(
    "-l",
    metavar="INT",
    type=int,
    dest="min_length",
    default=711,
    help="minimun length (with tag and primer) of COI barcode allowed, default=711",
)

polish_group.add_argument(
    "-L",
    metavar="INT",
    type=int,
    dest="max_length",
    default=719,
    help="maximun length (with tag and primer) of COI barcode allowed, default=719",
)

# ------------------------------------------------------------------------------------------------

###############################################################################
#####----------------------- main subcommand parsers --------------------######

description = """

Description

    An automatic pipeline for HIFI-SE400 project, including filtering
    raw reads, assigning reads to samples, assembly HIFI barcodes
    (COI sequences), polished assemblies, and do tax identification.
    See more: https://github.com/comery/HIFI-barcode-SE400

Versions

    1.0.5 (20190409)

Authors

    yangchentao at genomics.cn, BGI.
    mengguanliang at genomics.cn, BGI.
"""

parser = argparse.ArgumentParser(
    prog="HIFI-SE",
    description=description,
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument(
    "-v", "--version",
    action="version",
    version="%(prog)s 1.0.5"
)

subparsers = parser.add_subparsers(dest="command")

########## subcommmands ###########

## all subcommand
parser_all = subparsers.add_parser(
    "all",
    parents=[
        common_parser,
        index_parser,
        soft_parser,
        filter_parser,
        assign_parser,
        assembly_parser,
        trans_parser,
    ],
    formatter_class=argparse.RawTextHelpFormatter,
    help="run filter, assign and assembly.",
)

## filter subcommand
parser_filter = subparsers.add_parser(
    "filter",
    parents=[common_parser, filter_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="remove or trim reads with low quality.",
)

## assign subcommand
parser_assign = subparsers.add_parser(
    "assign",
    parents=[common_parser,
             index_parser,
             only_assign_parser,
             assign_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="assign reads to samples by tags.",
)

## assembly subcommand
parser_assembly = subparsers.add_parser(
    "assembly",
    parents=[common_parser,
             index_parser,
             only_assembly_parser,
             soft_parser,
             assembly_parser,
             trans_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="do assembly from assigned reads,\noutput raw HIFI barcodes.",
)

## polish subcommand
parser_polish = subparsers.add_parser(
    "polish",
    parents=[polish_parser,
            index_parser,
            trans_parser,],
    formatter_class=argparse.RawTextHelpFormatter,
    help="polish COI barcode assemblies,\n"
    + "output confident barcodes."
)

## BOLD_identification
parser_bold = subparsers.add_parser(
    "taxonomy",
    parents=[],
    formatter_class=argparse.RawTextHelpFormatter,
    help="do taxa identification on BOLD system\n",
)

###############################################################################
#####---------------------- program execution start ----------------------#####

# -----------------------BOLD identification----------------------#
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

if sys.argv[1] == "taxonomy":
    # if args.command == 'bold_identification':
    sys.argv = sys.argv[1:]
    sys.exit(bold_identification())

args = parser.parse_args()

# -----------------------arguments checking-----------------------#
## softwares and databases
def check_program_involed(cmd):
    '''
    check program involed whether is executable!
    '''
    result = (
        subprocess.call(
            "type %s" % cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        == 0
    )
    if result:
        return 0
    else:
        print(cmd + " not found!", file=sys.stderr)
        return 1


def files_exist_0_or_1(filelist):
    '''
    check files involed whether are existing!
    '''
    NUM = 0
    for file in filelist:
        if os.path.exists(file):
            NUM += 1
        else:
            print("%s doesn't exist!" % file, file=sys.stderr)
    if len(filelist) == NUM:
        return 0
    else:
        return 1

def print_time(str):
    print(str + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
# ----------------------------------------------------------------
## file existing check
errors_found = 0
if args.command == "all":
    errors_found += files_exist_0_or_1([args.raw, args.primer])
elif args.command == "filter":
    errors_found += files_exist_0_or_1([args.raw])
elif args.command == "assign":
    errors_found += files_exist_0_or_1([args.primer])
elif args.command == "assembly":
    errors_found += files_exist_0_or_1([args.list])
elif args.command == "polish":
    errors_found += files_exist_0_or_1([args.coi_input])
else:
    parser.print_help()
    parser.exit()

if args.command in ["all", "assembly"]:
    vsearch = "vsearch"
    if hasattr(args, "vsearch"):
        if args.vsearch:
            vsearch = args.vsearch
    errors_found += check_program_involed(vsearch)

if errors_found > 0:
    parser.exit("Errors found! Exit!")

if hasattr(args, "outpre") and args.outpre.endswith("/"):
    print("outpre is in bad format! no \"/\"")
    exit()

def check_and_open_outhandle(file):
    if os.path.exists(file):
        print("WARRNING: " + file + " exists! now overwriting ...")
    else:
        print("[INFO]: " + "open file " + file + "...")
    out = open(file, 'w')
    return out

# -----------------------functions for filtering------------------#



# ----------------------functions for assigning-------------------#

# ----------------------functions for assembling------------------#


# ------------------------filter process--------------------------#
# ------------------------assign process--------------------------#

if args.command in ["all", "assign"]:
    print_time("[INFO]: Assigning start:")


# ------------------------assembly process--------------------------#
if args.command in ["all", "assembly"]:
    print_time("[INFO]: Assembling start:")


#--------------------polish process-------------------------------------------
if args.command == "polish":
    print("[INFO]: Total barcodes polished: {}".format(polished_count))

print("[INFO]: total run time: {0:.2f}".format(time.time() - t) + "s")
