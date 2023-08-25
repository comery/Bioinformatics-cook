#!/usr/bin/env python3

import sys
import argparse


def addtwodict(thedict, ka, kb, val):
    if ka in thedict:
        thedict[ka].update({kb: val})
    else:
        thedict.update({ka: {kb: val}})

def check_exists(thedict, ka, kb):
    if ka in thedict:
        if kb in thedict[ka]:
            return True
        else:
            return False
    else:
        return False

def get_name(string):
    return string.split(".")[-2]


def main(args):
    seqid = set()
    info = {}
    score = {}
    databases = []
    if args.i:
        strings = args.i.strip().split(",")

    elif args.l:
        with open(args.l, 'r') as fh:
            strings = fh.readlines()
    else:
        sys.exit()
    out = open(args.o, 'w')
    for s in strings:
        database = get_name(s)
        databases.append(database)
        with open(s, 'r') as fh:
            for line in fh:
                tmp = line.strip().split("\t")
                seqid.add(tmp[0])
                addtwodict(info, tmp[0], database, "\t".join(tmp[1:]))
                addtwodict(score, tmp[0], database, float(tmp[11]))
    if args.m == 0:
        print(f"#SeqID\tDatabase\tInfo", file=out)
        for s in seqid:
            tmp_info = "-"
            tmp_score = 0
            for d in databases:
                if check_exists(score, s, d):
                    if score[s][d] > tmp_score:
                        tmp_score = score[s][d]
                        tmp_info = info[s][d]
            print(f"{s}\t[{d}]\t{tmp_info}", file=out)
    else:
        header = "\t".join(databases)
        print(f"#SeqID\t{header}", file=out)
        for s in seqid:
            tmp_out = []
            for d in databases:
                if check_exists(score, s, d):
                    tmp_out.append(info[s][d])
                else:
                    tmp_out.append("-")
            content = "\t".join(tmp_out)
            print(f"{s}\t{content}", file=out)
    out.close()

if __name__ == '__main__':
    des = """

    combine blast table into one table form several databases.
    each blast reuslt must have description on last column;
    for input file name, it needs to follow this format:
        xxx.xxx.nr.*; xxx.swissprot.tab; or something like that.

    yangchentao at genomics.cn, BGI.
    """
    parser = argparse.ArgumentParser(description=des,
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", type=str, required=False, metavar='<STR>',
                        help="input files, blast table format with dabase information")
    parser.add_argument("-l", type=str, required=False, metavar='<STR>',
                        help="input file list")
    parser.add_argument("-m", type=int, metavar='<INT>', default=1,
                        help="generating mode, 0 or 1, mixed or merged, default=1," +\
                       " mixed mode will output one item which have the highest score," +\
                       " merged mode will ouput all items in one line")
    parser.add_argument("-o", type=str, required=True, metavar='<STR>',
                        help="output file")
    args =  parser.parse_args()
    main(args)
    if not args.i and not args.l:
        sys.exit(f"you must specify one input")
