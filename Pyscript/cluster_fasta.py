import os
import argparse

parser = argparse.ArgumentParser(
	description="read fasta")

parser.add_argument('--fasta', metavar='FILE', type=str, required=True,
					help="input file, fasta file.")

args = parser.parse_args()

def match(str1,str2):
	base = 0
	matcched = 0
	while base <= len(str1):
		if str1[base] == str2[base]:
			matcched += 1
	return matcched / len(str1)


seq = {}

with open(args.fasta,'r') as fh:
	for line in fh:
		line = line.rstrip()
		if line[0] == '>':
			name=line[1:]
			seq[name] = ''
		else:
			seq[name]+= line.replace('\n','')

for k,v in seq.items():
	print(">" + k + "\n" + v )

sorted_name = sorted(seq.keys(), key=lambda x : len(seq[x]) , reverse = True )

tmp_seq = ""
for id in sorted_name:
	if tmp_seq != "":
				
