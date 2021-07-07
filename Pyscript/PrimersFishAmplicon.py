#!/usr/bin/python3
import os
import sys
import argparse
import subprocess
from Bio import SeqIO

usage = """
Description

    Extract amplicon sequence from your genome or long
    sequence by primer sets you given. you must give
    forward primer and reverse primer as 5'->3' order.


Usage

    python3 {0}  -f genome.fasta -p primers.fa -m 2

""".format(sys.argv[0])

parser = argparse.ArgumentParser(
    description=usage,
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-f', metavar='FILE', type=str, required=True,
                    dest='genome', help="your genome")

parser.add_argument('-p', metavar='FILE', type=str,
                   dest='primer', required=True,
                   help="primer fasta, formated as\n >TEST_F \n"
                   + "AGCTAGCTAGCTAGCTAG")

parser.add_argument('-m', metavar='INT', type=int,
                    dest='mismatch', default=2,
                    help="mismatch threshod for primer anchoring")

parser.add_argument('-o', metavar='STR', type=str,
                    dest='output', required=True,
                    help="output amplicons file")

if len(sys.argv) == 1:
    sys.exit(usage)

args = parser.parse_args()


# ------------------------------------------
def check_program_involed(cmd):
    '''
    check program involed whether is executable!
    '''
    result = (
        subprocess.call(
            "type %s" % cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        == 0
    )
    if result:
        return False
    else:
        return True

def parser_fasta(fh):
    name, seq = None, []
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def distance(s1, s2):
    code= {
            '-' : 0,
            'A' : 2,
            'T' : 3,
            'C' : 5,
            'G' : 7,
            'R' : 14,
            'Y' : 15,
            'M' : 20,
            'K' : 21,
            'S' : 35,
            'W' : 24,
            'H' : 30,
            'B' : 105,
            'V' : 70,
            'D' : 42,
            'N' : 210
    }

    match = 0
    for i in range(len(s1)):
        a = s1[i]
        b = s2[i]
        if a == b and a != "-" :
            match += 1
        elif code[a] + code[b] > 14 and code[b] % code[a] == 0:
            match += 1
        else:
            False
    mismatch = len(s1) - match
    return mismatch

def find_primer_binding(primer, seq, mismatch):
    potential_site = []
    plen = len(primer)
    qlen = len(seq)
    for i in range(qlen-plen):
        tmp = seq[i:i+plen].upper()
        if 'N' in tmp:
            continue
        else:
            dis = distance(tmp, primer)
            if dis <= mismatch:
                potential_site.append(i)
    if len(potential_site) >=2:
        print("more than 2 positions anchored!")
        #exit()
        return potential_site
    elif len(potential_site) == 0:
        return False
    else:
        return potential_site

def reverseComplement(sequence):
    # make a sequence complement #
    # replace function of string is too low!
    sequence = sequence.upper()
    transtable = str.maketrans('ATCG-', 'TAGC-')
    sequence = sequence.translate(transtable)
    return sequence

def addtwodimdict(thedict, key_a, key_b, val):
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})

#------------------------------------------
primers = {}
with open(args.primer, 'r') as fp:
    for i in parser_fasta(fp):
        name, seq = i
        name = name.replace(">", "")
        tags = name.split("_")
        if tags[1] == "F":
            addtwodimdict(primers, tags[0], tags[1], seq)
            revcom_F = reverseComplement(seq)
            addtwodimdict(primers, tags[0], 'Fc', revcom_F)
        else:
            addtwodimdict(primers, tags[0], tags[1], seq)
            revcom_R = reverseComplement(seq)
            addtwodimdict(primers, tags[0], 'Rc', revcom_R)

out = open(args.output, 'w')

with open(args.genome, 'r') as fg:
    for i in parser_fasta(fg):
        name, seq = i
        name = name.replace(">", "").split()[0]
        for pri in primers.keys():
            targetF = find_primer_binding(primers[pri]['F'], seq, args.mismatch)
            targetRc = find_primer_binding(primers[pri]['Rc'], seq, args.mismatch)
            targetR = find_primer_binding(primers[pri]['R'], seq, args.mismatch)
            targetFc = find_primer_binding(primers[pri]['Fc'], seq, args.mismatch)
            print(targetF, targetRc, targetR, targetFc)
            if targetF and targetRc:
                print("{0}, Primer {1}_F found".format(name, pri))
                print(targetF, targetRc)
                for x in targetF:
                    for y in targetRc:
                        if y > x:
                            length_Rc = len(primers[pri]['Rc'])
                            amplicon = seq[x:y + length_Rc]
                            out.write(">{0}_{1}_F_Rc;{2}_{3}\n{4}\n".format(name, pri, x, y, amplicon))
                        """
                        else:
                            length_F = len(primers[pri]['F'])
                            amplicon = seq[y:x + length_F]
                            out.write(">{0}_{1}_Rc_F;{2}_{3}\n{4}\n".format(name, pri, y, x, amplicon))
                        """

            elif targetR and targetFc:
                print("{0}, Primer {1}_R found".format(name, pri))
                print(targetR, targetFc)
                for x in targetR:
                    for y in targetFc:
                        if y > x:
                            length_Fc = len(primers[pri]['Fc'])
                            amplicon = seq[x:y + length_Fc]
                            amplicon = reverseComplement(amplicon)
                            out.write(">{0}_{1}_R_Fc;{2}_{3};reverseComplement\n{4}\n".format(name, pri, x, y, amplicon))
                        """
                        else:
                            length_R = len(primers[pri]['R'])
                            amplicon = seq[y:x + length_R]
                            out.write(">{0}_{1}_Fc_R;{2}_{3}\n{4}\n".format(name, pri, y, x, amplicon))
                        """

out.close()
