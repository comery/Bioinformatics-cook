#!/usr/env/python python3
import sys
import gzip
from icecream import ic

if len(sys.argv) != 4:
    print("Usage: python3 {} *.fastq1,*.fastq2 outprefix np_cutoff[int]".format(sys.argv[0]))
    exit()


def parse_se_fastq(fq_fh):
    """
    make generator to read fastq file.
    """
    while True:
        name = fq_fh.readline().strip()
        if len(name) == 0:
            break
        read = fq_fh.readline().strip()
        nothing = fq_fh.readline().strip()
        qual = fq_fh.readline().strip()
        yield name, read, qual


def smart_open(file):
    import gzip
    if file.endswith("gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file)

def main():
    np_cutoff = int(sys.argv[-1])
    outprefix = sys.argv[-2]
    low_pass_read = 0
    high_pass_read = 0
    low = gzip.open(outprefix + ".low.fasta.gz", 'wt')
    high = gzip.open(outprefix + ".high.fasta.gz", 'wt')
    stat = open(outprefix + ".filter.stat", 'w')
    for fq in sys.argv[1].split(","):
        with smart_open(fq) as fh:
            for i in parse_se_fastq(fh):
                name, read, qual = i
                # np=5 rq=0.995647
                name = name.replace("@", "")
                tmp = name.split()
                np = tmp[1].replace("np=", "")
                if int(np) < np_cutoff:
                    low.write(f">{name}\n{read}\n")
                    low_pass_read += 1
                else:
                    high.write(f">{name}\n{read}\n")
                    high_pass_read += 1
    total = low_pass_read + high_pass_read
    low_rate = "{:.3f}".format(low_pass_read/total)
    high_rate = "{:.3f}".format(high_pass_read/total)
    print(f"total ccs: {total}", file=stat)
    print(f"high quality ccs: {high_pass_read} ({high_rate})", file=stat)
    print(f"low quality ccs: {low_pass_read} ({low_rate})", file=stat)
    low.close()
    high.close()
    stat.close()

if __name__ == '__main__':
    main()


