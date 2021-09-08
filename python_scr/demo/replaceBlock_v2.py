#!/usr/bin/env python3
import sys
if len(sys.argv) < 4:
    sys.exit("python3 {} <*.1coords> <ref.fa> <querry.fa>".format(sys.argv[0]))


def get_id(seq_id):
    seq_id = seq_id.replace(">", "").split()[0]
    return seq_id

def parser_fasta(fasta):
    dict_fasta = {}
    with open(fasta, 'r') as fh:
        name, seq = None, []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    name = get_id(name)
                    dict_fasta[name] = "".join(seq)
                name, seq = line, []
            else:
                seq.append(line)
        if name:
            name = get_id(name)
            dict_fasta[name] = "".join(seq)
    return dict_fasta


def revcom(sequence):
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


def read_ref_block(inputfile):
    pre = ""
    with open(inputfile, 'r') as fh:
        while True:
            line = fh.readline().strip()
            if line.startswith("#"):
                continue
            if not line:
                break
            tmp = line.split()
            tmp = [tmp[0], int(tmp[1]), int(tmp[2]), tmp[3], int(tmp[4]), int(tmp[5]), tmp[-1]]
            if pre == "":
                block = [tmp,]
                pre = tmp[0]
            elif tmp[0] == pre:
                block.append(tmp)
            else:
                yield pre, block
                block = [tmp,]
                pre = tmp[0]

        yield tmp[0], block


def main():
    ref_seqs = parser_fasta(sys.argv[2])
    que_seqs = parser_fasta(sys.argv[3])
    replaced_scafs = []
    for rid, blocks in read_ref_block(sys.argv[1]):
        print(">" + rid + "_replaced")
        replaced_scafs.append(rid)
        pices = len(blocks)
        # replace whole scaffold
        if pices == 1 and blocks[0][-1] == "all":
            rid, rs, re, qid, qs, qe, rtype = blocks[0]
            if qs > qe:
                wrapseq(revcom(que_seqs[qid]))
            else:
                wrapseq(que_seqs[qid])
        # partially replace
        else:

            last_end = 0
            for b in blocks:
                rid, rs, re, qid, qs, qe, rtype = b
                # get aligned sequence
                # querry (hifi scaffold) strand = "+"
                if qe >= qs:
                    aligned = que_seqs[qid][qs-1:qe]
                # querry (hifi scaffold) strand = "-"
                else:
                    aligned = revcom(que_seqs[qid][qe-1:qs])

                # the first block
                if last_end == 0:
                    unaligned = ref_seqs[rid][0:rs-1] # unaligned sequence between two blocks of reference
                    makeup_superscaf  = unaligned + aligned
                else:
                    # overlapped
                    if rs <= last_end:
                        overlapped_length = last_end - rs + 1
                        # trim the end of makeup_superscaf
                        makeup_superscaf = makeup_superscaf[:-overlapped_length]
                        makeup_superscaf += aligned
                    else:
                        unaligned = ref_seqs[rid][last_end:rs-1]
                        makeup_superscaf  += unaligned + aligned
                last_end = re

            # deal with tail
            if last_end < len(ref_seqs[rid]):
                makeup_superscaf += ref_seqs[rid][last_end:]
            wrapseq(makeup_superscaf)

    # output other scaffolds

    for k in ref_seqs.keys():
        if k not in replaced_scafs:
            print(">" + k)
            wrapseq(ref_seqs[k])


if __name__ == '__main__':
    main()
