import sys
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def wrap(seq):
    tmp = 0
    length= len(seq)
    new = ""
    while tmp < length:
        new += seq[tmp:tmp+100] + "\n"
        tmp += 100
    return new

with open(sys.argv[1] ,'r') as fp:
    for name, seq in read_fasta(fp):
        print(name, seq)
