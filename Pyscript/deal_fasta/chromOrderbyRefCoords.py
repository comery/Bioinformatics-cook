#!/usr/bin/env python3
import sys
if len(sys.argv) < 4:
    sys.exit("python3 {} <*.1coords> <ref: scaf2chr.bed> <querry.fa>".format(sys.argv[0]))


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
    step = 100
    new = ""
    for i in range(0,l,step):
        if i+step <= l:
            new += seq[i:i+step] + "\n"
        else:
            new += seq[i:]

    return new


def addtwodict(thedict, k1, k2, val):
    if k1 in thedict:
        thedict[k1].update({k2: val})
    else:
        thedict.update({k1: {k2: val}})


def read_coords(inputfile):
    coords = {}
    que_lens = {}
    with open(inputfile, 'r') as fh:
        while True:
            line = fh.readline().strip()
            if line.startswith("#"):
                continue
            if not line:
                break
            items = line.split()
            # convert str type to int
            start1,end1,start2,end2,len1,len2,identity,len_chr1,len_chr2,cov_chr1,cov_chr2,chr1,chr2 = items
            [start1,end1,start2,end2,len1,len2,len_chr1,len_chr2] = [ int(i) for i in [start1,end1,start2,end2,len1,len2,len_chr1,len_chr2] ]
            [identity, cov_chr1, cov_chr2] = [ float(f) for f in [identity, cov_chr1, cov_chr2] ]
            # query strand
            if start2 < end2:
                strand = 1
            else:
                strand = -1
            # store each coord relation to 2 dim dict, key1 is ref_scaf_id, key2 is start position, val is a dict containing all infomation
            tmp = {'chr1':chr1,'start1':start1,'end1':end1,'chr2':chr2,'start2':start2,'end2':end2,'cov_chr1':cov_chr1,'cov_chr2':cov_chr2,'strand':strand}
            addtwodict(coords, chr1, start1, tmp)
            # store ref scaffold length for sorting them by length later
            que_lens[chr2] = len_chr2
    return coords, que_lens


class block:
    def __init__(self, record):
        self.chr1 = record['chr1']
        self.start1 = record['start1']
        self.end1 = record['end1']
        self.cov1 = record['cov_chr1']
        self.chr2 = record['chr2']
        self.start2 = record['start2']
        self.end2 = record['end2']
        self.cov2 = record['cov_chr2']
        self.strand = record['strand']
    def update(self, key, val):
        if key not in ['start1', 'end1', 'start2', 'end2', 'strand']:
            sys.exit("can not accept this key to update")
        else:
            self.key = val


def get_related_blocks(tmp_dict):
    all_pos = sorted(tmp_dict.keys())
    pre = ""
    related_blocks = []
    for p in all_pos:
        b = block(tmp_dict[p])
        if pre == "":
            pre = b.chr2
            related_blocks.append(b)
        elif b.chr2 == pre:
            related_blocks.append(b)
        else:
            yield related_blocks
            related_blocks = [b,]
            pre = b.chr2

    yield related_blocks

def merge_blocks(related_blocks):
    sum_strand = 0
    chr2_poses = []
    chr1_poses = []
    for b in related_blocks:
        chr2_poses.append(b.start2)
        chr2_poses.append(b.end2)
        chr1_poses.append(b.start1)
        chr1_poses.append(b.end1)
        sum_strand += b.strand
    chr2_poses.sort()
    chr1_poses.sort()
    b.update('start1', chr1_poses[0])
    b.update('end1', chr1_poses[-1])
    b.update('start2', chr2_poses[0])
    b.update('end2', chr2_poses[-1])
    b.update('strand', sum_strand)

    return b


def main():

    que_seqs = parser_fasta(sys.argv[3])
    replaced_scafs = []
    new_que_seq = []
    que_bed = []
    current_que_start = 1
    coords, que_lens = read_coords(sys.argv[1])

    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            if i.startswith("#"):continue
            tmp = i.strip().split()
            if tmp[-2] == 'gap':continue

            ref_scaf = tmp[-2]
            [ref_scaf_start, ref_scaf_end] = [int(tmp[1]), int(tmp[2])]
            ref_scaf_len = ref_scaf_end - ref_scaf_start + 1

            if ref_scaf in coords:
                # related blocks are from the same query and they are linked in location
                for related_blocks in get_related_blocks(coords[ref_scaf]):
                    # only one block
                    if len(related_blocks) < 2:
                        b = related_blocks[0]
                    else:
                        # merge related blocks into one block, and decide the strand
                        b = merge_blocks(related_blocks)

                    if b.strand >= 0:
                        new_que_seq.append(que_seqs[b.chr2])
                        current_que_end = current_que_start + que_lens[b.chr2] - 1
                        que_bed.append(f"Que_chr1\t{current_que_start}\t{current_que_end}\t{b.chr2}\t+")
                    else:
                        new_que_seq.append(revcom(que_seqs[b.chr2]))
                        current_que_end = current_que_start + que_lens[b.chr2] - 1
                        que_bed.append(f"Que_chr1\t{current_que_start}\t{current_que_end}\t{b.chr2}\t-")
                    current_que_start = current_que_end + 1

                new_que_seq.append("n" * 500)
                current_que_end = current_que_start + 500 - 1 # add gap
                que_bed.append(f"Que_chr1\t{current_que_start}\t{current_que_end}\tglue\t+")
                current_que_start = current_que_end + 1
            else:
                if len(new_que_seq) > 0:
                    new_que_seq.pop()
                    que_bed.pop()
                new_que_seq.append("n" * ref_scaf_len)
                current_que_end = current_que_start + ref_scaf_len - 1 # add gap
                que_bed.append(f"Que_chr1\t{current_que_start}\t{current_que_end}\tgap:{ref_scaf}\t+")
                current_que_start = current_que_end + 1


    with open("new_que.fa", 'w') as fh:
        print(">Que_chr1", file=fh)
        print(wrapseq("".join(new_que_seq)), file=fh)

    with open("new_que.bed", 'w') as bh:
        print("\n".join(que_bed), file=bh)



if __name__ == '__main__':
    main()
