#!/usr/bin/python3
import sys
import argparse

def revcom(sequence):
    # make a sequence complement #
    sequence = sequence.upper()
    transtable = str.maketrans('ATCG', 'TAGC')
    sequence = sequence.translate(transtable)
    return sequence[::-1]

def search(fa, seq, tag):
    if fa.endswith('.gz'):
        import gzip
        f = gzip.open(fa)
    else:
        f = open(fa)
    import re
    from Bio import SeqIO
    out = open("telomere." + tag + ".bed", 'w')
    pattern = re.compile(seq)
    for record in SeqIO.parse(f, 'fasta'):
        for m in pattern.finditer(str(record.seq).lower()):
            print(f"{record.id}\t{m.start()+1}\t{m.end()}", file=out)
    out.close()


def main():

    parser = argparse.ArgumentParser(description='genomic telomere scan tool')
    parser.add_argument('faFile', type=str, help="input fasta")
    parser.add_argument('motif', type=str, help='telomere motif [e.g. TTAGGG]')
    parser.add_argument('strand', type=int, choices=[0,1,2], help='search type, 0:itself, 1:reverse & complement, 2:both')

    args = parser.parse_args()

    faFile = args.faFile
    seq = args.motif.lower()
    seq1 = revcom(seq).lower()
    #print(f"{seq}\t{seq1}")


    if args.strand == 0 or args.strand == 2:
        search(faFile, seq, '0')
    if args.strand == 1 or args.strand == 2:
        search(faFile, seq1, '1')

if __name__ == '__main__':
    main()

