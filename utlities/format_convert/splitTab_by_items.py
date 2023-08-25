#!/usr/bin/env python3
import sys
import os

if len(sys.argv) < 4:
    sys.exit(f"python3 {sys.argv[0]} <tsv> <column,1-based> <outdir>")

def get_key(line, column):
    items = column.split(",")
    keys = []
    for k in items:
        try:
            k = int(k)
        except TypeError:
            print("please give integer!")

        if k < 1:
            sys.exit("can be 0 for column!")
        keys.append(line[k-1])
    return "-".join(keys)


def main():
    records = {}
    outdir = sys.argv[3]
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            i = i.strip()
            if i.startswith("#"):continue
            tmp = i.strip().split()
            key = get_key(tmp, sys.argv[2])
            if key in records:
                records[key].append(i)
            else:
                records[key] = [i,]


    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    for k in records:
        outfile = os.path.join(outdir, k + ".txt")
        with open(outfile, 'w') as fh:
            print("\n".join(records[k]), file=fh)

    print("All done!")

if __name__ == '__main__':
    main()
