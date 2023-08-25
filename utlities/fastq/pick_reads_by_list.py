#!/usr/env/python python3
import sys

if len(sys.argv) != 4:
    print("Usage: python3 {} *.fastq1 *.fastq2 *.list".format(sys.argv[0]))
    exit()


def parse_se_fastq(fq_fh):
    """
    make generator to read fastq file.
    """
    while True:
        names = fq_fh.readline().strip().split()
        if len(names) == 0:
            break
        name = names[0]
        name = name[1:]
        read = fq_fh.readline().strip()
        nothing = fq_fh.readline().strip()
        qual = fq_fh.readline().strip()
        yield name, read, qual


list = []
with open(sys.argv[3], 'r') as fh:
    for i in fh:
        list.append(i.strip())

#print(list)


clean_fq1 = open("clean_1.mito.fq", 'w')
clean_fq2 = open("clean_2.mito.fq", 'w')

with open(sys.argv[1], 'r') as fh:
   for i in parse_se_fastq(fh):
        name, read, qual = i
        #print(name)
        if name in list:
            clean_fq1.write("@" + name + "\n" + read + "\n" + "+\n" + qual + "\n")

#clean_fq1.close()


with open(sys.argv[2], 'r') as fh:
    for i in parse_se_fastq(fh):
        name, read, qual = i
        #print(name)
        if name in list:
            clean_fq2.write("@" + name + "\n" + read + "\n" + "+\n" + qual + "\n")

clean_fq2.close()

