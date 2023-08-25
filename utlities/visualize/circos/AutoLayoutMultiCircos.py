#!/usr/bin/env python3

import sys

if len(sys.argv) < 4:
    sys.exit(f"python3 {sys.argv[0]} <min_radious> <max_radious> <list>")

def main():
    min_radious = float(sys.argv[1])
    max_radious = float(sys.argv[2])
    files = []
    with open(sys.argv[3], 'r') as fh:
        for i in fh:
            if i.startswith("#"):continue
            files.append(i.strip())
    num = len(files)
    r = 5
    aa = (max_radious - min_radious)/(r*num + num - 1)
    bb = r * aa
    unit = aa + bb;

    files.sort()
    for i, f in enumerate(files):
        r0 = min_radious + i * unit
        r1 = r0 + bb
        print("<highlight>")
        print(f"file    = {f}")
        print(f"r0      = {r0}r")
        print(f"r1      = {r1}r")
        print(f"fill_color  = gray")
        print("</highlight>")
    print("</highlights>")

if __name__ == '__main__':
    main()
