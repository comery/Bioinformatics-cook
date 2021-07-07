#!/usr/bin/env python3
import sys
if len(sys.argv) < 3:
    sys.exit("python3 {} <scaf.fa> <*.agp>".format(sys.argv[0]))


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, "".join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, "".join(seq))


def revcom_transtable(sequence):
    # make a sequence complement #
    transtable = str.maketrans('ATCGatcg', 'TAGCtagc')
    sequence = sequence.translate(transtable)
    return sequence[::-1]

def wrapseq(seq):
    l = len(seq)
    for i in range(0,l,80):
        if i+80 <= l:
            print(seq[i:i+80])
        else:
            print(seq[i:])

def main():
    scafs = {}
    with open(sys.argv[1], 'r') as fa:
        for name, seq in read_fasta(fa):
            name = name.replace(">", "")
            scafs[name] = seq
    with open(sys.argv[2], 'r') as fh:
        data = {}
        superscafs = set()
        for i in fh:
            tmp = i.strip().split()
            if len(tmp) != 9:
                sys.exit("bad format of {}".format(sys.argv[2]))
            superscafs.add(tmp[0])
            if tmp[0] in data.keys():
                data[tmp[0]].append(tmp)
            else:
                data[tmp[0]] = [tmp,]

    for s in list(superscafs):
        makeup_superscaf = ""
        for r in data[s]:
            pos_s = int(r[1])
            pos_e = int(r[2])
            len1 = pos_e - pos_s + 1
            comp_type = r[4]
            subscaf = r[5]
            subscaf_beg = int(r[6])
            subscaf_end = int(r[7])
            len2 = subscaf_end - subscaf_beg + 1
            strand = r[8]
            if len1 != len2:
                sys.exit("wrong position of {}".format(sys.argv[2]))

            if strand == "+":
                if comp_type == 'W':
                    makeup_superscaf += scafs[subscaf][subscaf_beg-1:subscaf_end]
                elif comp_type == 'N':
                    makeup_superscaf += 'N'*subscaf_end
                else:
                    print("this type {} is not accepted".format(comp_type))
            else:
                if comp_type == 'W':
                    makeup_superscaf += revcom_transtable(scafs[subscaf][subscaf_beg-1:subscaf_end])
                elif comp_type == 'N':
                    makeup_superscaf += 'N'*subscaf_end
                else:
                    print("this type {} is not accepted".format(comp_type))
        print(">" + s)
        wrapseq(makeup_superscaf)

if __name__ == '__main__':
    main()
