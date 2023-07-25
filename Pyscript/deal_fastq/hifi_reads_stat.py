#!/usr/env/python python3
import sys
import pyfastx

if len(sys.argv) != 2:
    print("Usage: python3 {} *.fastq1,*.fastq2".format(sys.argv[0]))
    exit()


def store_to_db(data, outfile):
    with open(outfile, 'w') as fw:
        fw.write(json.dumps(data, indent=4))

def statistic(filename, lens):
    print(f"# {filename}")
    # summary and stat
    total_bases = sum(lens)
    total_ccs = len(lens)
    print(f"total bases = {total_bases}")
    print(f"total ccs = {total_ccs}")


def main():
    read_bases = []
    for fa in sys.argv[1].split(","):
        current_file = []
        for seq in pyfastx.Fasta(fa):
            current_file.append(len(seq.seq))
            read_bases.append(len(seq.seq))
        statistic(fa, current_file)

    # overall
    statistic("Above All", read_bases)

if __name__ == '__main__':
    main()


