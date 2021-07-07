import sys
try:
    import Bio
except:
    sys.exit("package biopython not found! Please install it!")
else:
    #from Bio.Alphabet import IUPAC
    from Bio import Seq

def translate_dnaseq(seq, codon):
# ---------translate_dnaseq------------#
    l_dna = len(seq)
    if l_dna % 3 != 0:
        print("DNA sequence length is not triplet!")
    coding_dna = Seq.Seq(seq)
    protein = coding_dna.translate(table=codon)
    return protein


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


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(print(f"python3 {sys.argv[0]} <fasta> <codon table>"))

    codon = sys.argv[2]
    with open(sys.argv[1], 'r') as fp:
        for name, seq in read_fasta(fp):
            pro = translate_dnaseq(seq, codon)
            print(name + "\n" + pro)
