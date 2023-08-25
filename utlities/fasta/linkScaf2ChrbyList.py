#!/usr/bin/env python3
import sys
if len(sys.argv) < 4:
    sys.exit("python3 {} <assembly.fa> <scaffold list> <Num of N to glue>".format(sys.argv[0]))


def get_id(seq_id):
    seq_id = seq_id.replace(">", "").split()[0]
    return seq_id

def parser_fasta(fasta):
    seqs = {}
    scaf_lens = {}
    with open(fasta, 'r') as fh:
        name, seq = None, []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    name = get_id(name)
                    seqs[name] = "".join(seq)
                    scaf_lens[name] = len("".join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name:
            name = get_id(name)
            seqs[name] = "".join(seq)
            scaf_lens[name] = len("".join(seq))
    return seqs, scaf_lens


def wrapseq(seq):
    l = len(seq)
    step = 100
    new = ""
    for i in range(0,l,step):
        if i+step <= l:
            new += seq[i:i+step] + "\n"
        else:
            new += seq[i:]

    return new

def main():
    seqs,scaf_lens = parser_fasta(sys.argv[1])
    Nglue = int(sys.argv[3])
    glue = "n" * Nglue
    replaced_scafs = []
    new_ref_seq = []
    out1 = open("scaf2chr.bed", 'w')
    current_ref_start = 1
    scaf_order = []
    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            if i.startswith("#"):continue
            tmp = i.strip()
            scaf_order.append(tmp)
    for ref_chr in scaf_order:
        # store each ref sequence to a list for output
        new_ref_seq.append(seqs[ref_chr])
        # to generate agp file for newly makeuped ref seq
        current_ref_end = current_ref_start + scaf_lens[ref_chr] - 1
        print(f"Ref_chr1\t{current_ref_start}\t{current_ref_end}\t{ref_chr}\t+", file=out1)
        current_ref_start = current_ref_end + 1
        # related blocks are from the same query and they are linked in location

        if ref_chr != scaf_order[-1]: #don't add gap at the end
            current_ref_end = current_ref_start + Nglue - 1 # add gap
            print(f"Ref_chr1\t{current_ref_start}\t{current_ref_end}\tgap\t+", file=out1)
            current_ref_start = current_ref_end + 1

    out1.close()

    print(">Ref_chr1")
    print(wrapseq(glue.join(new_ref_seq)))


if __name__ == '__main__':
    main()
