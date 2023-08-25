import sys
if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} *.pep")
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, "\n".join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, "\n".join(seq))

with open(sys.argv[1], 'r') as fp:
    for name, seq in read_fasta(fp):
        if 'U' not in seq:
            print(name + "\n" + seq)
