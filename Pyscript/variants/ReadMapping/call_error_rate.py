#!/usr/bin/env python3
import sys
import pyfastx
from icecream import ic


def parse_fasta(faFile):
    dic = {}
    fa = pyfastx.Fastx(faFile)
    for name,seq,comment in fa:
        dic[name] = seq.strip().replace("\n", "")

    return dic


def main():
    ref = parse_fasta(sys.argv[2])
    callable_site = 0
    ins_site = 0
    del_site = 0
    error_site = 0
    callable_base = 0
    error_base = 0
    genotype_ins = {}
    genotype_del = {}
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i.startswith("Scaffold"):
                continue
            tmp = i.strip().split("\t")
            scaf = tmp[0]
            loca = int(tmp[1])
            depth = int(tmp[2])
            a, t, g, c, pi, ins, dels = tmp[4:11]
            a = int(a)
            t = int(t)
            g = int(g)
            c = int(c)
            pi = float(pi)
            if depth < 3:
                continue
            callable_site += 1
            # callable base
            callable_base += depth
            tmp_dic = {'A': a,
                       'T': t,
                       'G': g,
                       'C': c}
            ref_base = ref[scaf][loca-1]
            if ref_base == "N":
                continue
            #ic(scaf, loca, ref_base)
            del tmp_dic[ref_base]
            error_base += sum(tmp_dic.values())

            if pi > 0 :
                error_site += 1
            elif ins:
                insertions = ins.split(";")
                for x in insertions:
                    genoty = x.split("(")[0]
                    if "N" in genoty:
                        continue
                    ins_site += 1
                    if genoty in genotype_ins:
                        genotype_ins[genoty] += 1
                    else:
                        genotype_ins[genoty] = 1
            elif dels:
                deletions = dels.split(";")
                for y in deletions:
                    genoty = y.split("(")[0]
                    if "N" in genoty:
                        continue
                    del_site += 1
                    if genoty in genotype_del:
                        genotype_del[genoty] += 1
                    else:
                        genotype_del[genoty] = 1
    error_rate_site  = "{:.5f}".format(error_site/callable_site)
    error_rate_base  = "{:.5f}".format(error_base/callable_base)
    accuracy = "{:.5f}".format(1-error_base/callable_base)
    print(f"callable site: {callable_site}")
    print(f"error site: {error_site}")
    print(f"error rate (site): {error_rate_site}")
    print(f"callable base: {callable_base}")
    print(f"error base: {error_base}")
    print(f"accuracy: {accuracy}")
    print(f"error rate (base): {error_rate_base}")
    print(f"insertion site: {ins_site}")
    print(f"deletion site: {del_site}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(f"python3 {sys.argv[0]} base.xls ref")
    else:
        main()
