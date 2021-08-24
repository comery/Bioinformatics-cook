#!/usr/bin/env python3
import sys
import re

if len(sys.argv) < 3:
    sys.exit("python3 {} <codon.bed> <sites>".format(sys.argv[0]))


def main():
    psgs = {}
    genemap = {}
    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            tmp = i.strip().split("\t")
            genename = tmp[0]
            pid = tmp[1]
            genemap[pid] = genename
            if pid in psgs.keys():
                psgs[pid].append(tmp[-2])
            else:
                psgs[pid] = [tmp[-2],]

# MARMOS_MARMOS_R00001  Super_scaffold_mat_1    +   84904361    84904362    84904363    1
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            pid = tmp[0].replace("MARMOS_MARMOS", "MARMOS")
            codon = tmp[-1]
            if pid in psgs.keys():
                if codon in psgs[pid]:
                    print("{}\t{}\t{}\t{}".format(genemap[pid], pid, codon, i.strip()))


if __name__ == '__main__':
    main()
