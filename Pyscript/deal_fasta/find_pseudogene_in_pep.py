import sys
if len(sys.argv) < 3:
    sys.exit("Usage: python3 {} <input.pep> <out.pep>".format(sys.argv[0]))

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


pep = {}
marm_stop = {}
out = open(sys.argv[2], 'w')
pseudogene = 0
with open(sys.argv[1], 'r') as fp:
    for name, seq in read_fasta(fp):
        name = name.replace(">", "")
        if 'U' in seq:
            pseudogene = 1
        else:
            print(">" + name + "\n" + seq, file=out)

out.close()
if pseudogene > 0:
    print(sys.argv[1] + " pseudogene")


