#!/usr/bin/env python3
import sys
if len(sys.argv) < 2:
    print("python3 {} <*.txt>".format(sys.argv[0]))
    exit()

import re
rgb_re = re.compile(r'^(\d+),(\d+),(\d+)$')
with open(sys.argv[1]) as fh:
    for i in fh:
        items = i.strip().split()
        for it in items:
            new = []
            m = rgb_re.match(it)
            if m:
                tmp = it.split(",")
                r, g, b = int(tmp[0]), int(tmp[1]), int(tmp[2])
                new.append('#{:02x}{:02x}{:02x}'.format(r, g, b))
            else:
                new.append(it)

            print("\t".join(new))

