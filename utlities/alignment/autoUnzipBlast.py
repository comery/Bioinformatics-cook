#!/usr/bin/env python3
import os
import sys
import subprocess
import argsparser

def open_input(file):
    if file.endswith("gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')

"A R N D C E Q G H I L K M F P S T W Y V"
def detectFastaType(fasta):
    with open(fasta, 'r') as fh:
        name = fh.readline().strip()
        seq = fh.readline().strip()[:100]
        if not name.startswith(">"):
            print("This is not a FASTA file!")
            exit()
        elif 




def makeblastdb()
