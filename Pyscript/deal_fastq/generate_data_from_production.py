#/usr/bin/env python3
import sys
import os
import time
import glob
import argparse
from icecream import ic


def find_max_order(outdir, symb):
    found = glob.glob(outdir + symb)
    if len(found) == 0:
        return 0
    else:
        print(f"find existing files {found}")
        return len(found)


def main():
    date = time.strftime("%b-%d", time.localtime())
    records = {}
    rawpaths = {}
    sample2individual = {}
    with open(args.map, 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            sample2individual[tmp[0]] = tmp[1]

    # get a non-redundant list of raw paths
    with open(args.input, 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            raw_path, sample = tmp[4], tmp[-3]
            if sample not in rawpaths:
                rawpaths[sample] = set()
                rawpaths[sample].add(raw_path)
            else:
                rawpaths[sample].add(raw_path)
    # get a non-redundant list of (fq1, fq2)
    for sample in rawpaths.keys():
        for tmp_path in list(rawpaths[sample]):
            for fq in glob.glob(tmp_path + "/*_1.fq.gz"):
                fq2 = fq.replace("_1.fq.gz", "_2.fq.gz")
                abs_fq1 = os.path.join(tmp_path, fq)
                abs_fq2 = os.path.join(tmp_path, fq2)
                if sample not in records:
                    records[sample] = [(abs_fq1, abs_fq2),]
                else:
                    records[sample].append((abs_fq1, abs_fq2))
    if args.filter:
        target_path = args.filter
        if os.path.exists(target_path) == False:
            sys.exit(f"no such dir: {target_path}")
        for sample in records:
            indi =  sample2individual[sample]
            print(f"dealing with filter of {indi}...")
            outdir = os.path.join(os.path.abspath(target_path), indi)
            if os.path.exists(outdir) == False:
                os.mkdir(outdir)
            existing_max_order = find_max_order(outdir, "/*_1.clean.fq.gz")
            filter_shell = os.path.join(outdir, "filter" + "_" + date + ".sh")
            with open(filter_shell, 'w') as oh:
                for index, fq in enumerate(records[sample]):
                    order = index + existing_max_order + 1
                    fq1, fq2 = fq
                    out1 = os.path.join(outdir, indi + "-" + str(order) + "_1.clean.fq.gz")
                    out2 = os.path.join(outdir, indi + "-" + str(order) + "_2.clean.fq.gz")
                    out_json = os.path.join(outdir, indi + "-" + str(order) + ".json")
                    out_html = os.path.join(outdir, indi + "-" + str(order) + ".html")
                    print("date", file=oh)
                    print(f"fastp -i {fq1} -o {out1} -I {fq2} -O {out2} -l 60 -u 20 -w {args.threads} -j {out_json} -h {out_html}", file=oh)
                    print("date", file=oh)
    # backup
    if args.backup:
        #ic(args.backup)
        backup_path = args.backup
        if os.path.exists(backup_path) == False:
            sys.exit(f"can not find backup path: {backup_path}")
        for sample in records:
            indi =  sample2individual[sample]
            print(f"dealing with backup of {indi}...")
            outdir = os.path.join(os.path.abspath(backup_path), indi)
            if os.path.exists(outdir) == False:
                os.mkdir(outdir)
            backup_shell = os.path.join(outdir, "backup" + "_" + date + ".sh")
            #ic(backup_shell)
            with open(backup_shell, 'w') as bh:
                for index, fq in enumerate(records[sample]):
                    fq1, fq2 = fq
                    raw_dir = os.path.dirname(fq1)
                    print("date", file=bh)
                    print(f"cp {raw_dir}/*  {outdir}/ && echo 'done' ", file=bh)
                    print("date", file=bh)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='input an og draw the positively selected sites for manual check of reliable sites and alignment issue sites')
    parser.add_argument('input', type=str, help="data.list")
    parser.add_argument('-b', '--backup', type=str, dest='backup', help='a writable path for backup')
    parser.add_argument('-m', '--map', type=str, dest='map', help='tab that contains sample and corresponding individual')
    parser.add_argument('-f', '--filter', type=str, dest='filter', help='a writable path for clean data')
    parser.add_argument('-t', '--threads', type=int, default=3, dest='threads', help='number of threads for fastp, default=3')
    parser.add_argument('-s', '--stat', dest='stat', action="store_true", help="get statment of data")

    args = parser.parse_args()
    main()
