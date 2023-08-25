import sys
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


if len(sys.argv) < 2:
    print("python3 {} fasta".format(sys.argv[0]))
    exit()
with open(sys.argv[1], 'r') as fh:
    for i in read_fasta(fh):
        name, seq = i
        print(f"{name}\t{len(seq)}")


