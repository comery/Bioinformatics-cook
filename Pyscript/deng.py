import os
import sys

usage = "usage:\n\tpython3 " + sys.argv[0] + " file1,file2,file3... " + " need_search file"

print(len(sys.argv))
if len(sys.argv) < 3:
    print(usage)
    exit()

o_list = sys.argv[1]
s_file = sys.argv[2]

pairs = []

o_lists = o_list.split(',')
print(o_lists)

for f in o_lists:
	with opne(f) as fh_f:
		obj = fh_f.readlines()
		for line in obj:
			gene = line.split('\t',2)
			print(gene)
			pairs.append(gene)

print("reading files done!" +\
	"now searching genes in " + s_file)


fh_out = open("output.txt", 'w')

with open(s_file) as fh_s:
	for line in fh_s.readline():
		tag = line.split('\t',2)
		if tag in pairs:
			print(line)

