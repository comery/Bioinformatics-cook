#!/usr/bin/env python3
import os
import sys
"""
This is for dengyuan, MHC analysis.
"""
if len(sys.argv) < 3:
    help="python3 {0} aligned.fasta coverage[INT]".format(sys.argv[0])
    print(help)
    exit()

if sys.argv[2].isdigit() == False:
    print("coverage is a int")
    exit()
else:
    deep_coverage = int(sys.argv[2])


seqs = {}
with open(sys.argv[1], 'r') as fh:
    for line in fh:
        if line.startswith('>'):
            line = line.strip()
            name=line.replace('>', '').split()[0]
            if name in seqs.keys():
                print("Error: " + name + " appears more than once!")
                exit()
            else:
                seqs[name] = ''

        else:
            seqs[name]+= line.replace('\n','')

max_id_len = 0
max_seq_len = 0
for k in sorted(seqs.keys()):
    if len(k) > max_id_len:
        max_id_len = len(k)
    if max_seq_len == 0:
        max_seq_len = len(seqs[k])
    else:
        if len(seqs[k]) != max_seq_len:
            print("Error: your file " + sys.argv[1] + " is not aligned!")
            exit()

with open(sys.argv[1] + ".link", 'w') as oh:
    for k in sorted(seqs.keys()):
        fill = " " * (max_id_len - len(k) + 5)
        oh.write(k + fill + seqs[k] + "\n")

records = len(seqs.keys())
conserved_blocks = []
start = 0
end = 0

for s in range(max_seq_len):
    tmp = [n[s] for n in seqs.values()]
    N_count = tmp.count('-')
    if N_count < records * 0.2:
        if s == 0:
            start = 1
            end = 1
        elif start == 0:
            start = s
            end = s
        elif start != 0:
            end += 1
    else:
        if start == 0:
            continue
        elif end - start + 1 >= deep_coverage:
            conserved_blocks.append((start, end))
            start = 0
        else:
            start = 0
            continue

with open(sys.argv[1] + ".blocks", 'w') as fh:
    for k in seqs.keys():
        fh.write(">" + k + " ")
        blocks = ""
        rang = ""
        for i in conserved_blocks:
            (a, b) = i
            t = "[{0}, {1}]".format(a, b)
            rang += t + " "
            block = seqs[k][a:b+1]
            blocks += block + " "
        fh.write(rang + "\n")
        fh.write(blocks + "\n")

