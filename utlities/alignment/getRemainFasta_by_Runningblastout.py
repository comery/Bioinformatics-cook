#!/usr/bin/env python3

import sys

if len(sys.argv) < 4:
    print("Usage: python3 {} {} {} {}".format(sys.argv[0], "<input.fasta>", "<stoped ID>", "<output>"))
    exit()

stop_id = sys.argv[2]


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name.replace(">", ""), '\n'.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name.replace(">", ""), '\n'.join(seq))


with open(sys.argv[1], 'r') as fh, open(sys.argv[3], 'w') as fw:
    sign = 0
    for name,seq in read_fasta(fh):
        ids = name.split()[0]
        #print(name + "\t" + ids)
        if ids == stop_id:
            print(">{}\n{}".format(name, seq), file=fw)
            sign = 1

        if sign == 0:
            continue
        else:
            print(">{}\n{}".format(name, seq), file=fw)


