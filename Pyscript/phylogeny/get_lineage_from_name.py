#!/usr/bin/python3
import sys
from ete3 import NCBITaxa
import argparse

description = """
output lineage of input taxonomy name.
"""
parser = argparse.ArgumentParser(description=description)

parser.add_argument("-in", dest="taxa_file", metavar="<file>",
					help="input taxonomy name list. One name" +\
                    "(e.g. species name) per-line, words" +\
                    "of one name should be seperated by a blank.")

parser.add_argument("-out", dest="out_file", metavar="<file>", help="outfile")

parser.add_argument("-lineage", dest="lineage_file", metavar="<file>",
					help="also output lineage tree (default: No)")

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit()

else:
	args = parser.parse_args()

name_list = []
with open(args.taxa_file, 'r') as fh:
	for i in fh:
		i = i.strip()
		if "_" in i:
			i = i.replace("_", " ")
		name_list.append(i)

ncbi = NCBITaxa()
taxi_name_dict = ncbi.get_name_translator(name_list)

taxi_id_list = [taxi_name_dict[i][0] for i in taxi_name_dict.keys()]

print(len(taxi_id_list), "of input taxa found.")

fh_out= open(args.out_file, 'w')

for taxi_id in taxi_id_list:
	taxi_id_lineage = ncbi.get_lineage(taxi_id)

	rank_dict = dict()
	for rank in ['kingdom', 'phylum', 'order', 'family', 'genus', 'species']:
		rank_dict[rank] = 'NA'

	for j in taxi_id_lineage:
		rank = ncbi.get_rank([j])[j]
		taxa = ncbi.get_taxid_translator([j])[j]
		if rank == 'kingdom':
			rank_dict['kingdom'] = taxa

		elif rank == 'phylum':
			rank_dict['phylum'] = taxa

		elif rank == 'order':
			rank_dict['order'] = taxa

		elif rank == 'family':
			rank_dict['family'] = taxa

		elif rank == 'genus':
			rank_dict['genus'] = taxa

		elif rank == 'species':
			rank_dict['species'] = taxa

		else:
			pass

	for rank in ['kingdom', 'phylum', 'order', 'family', 'genus', 'species']:
		print(rank_dict[rank], end="|", file=fh_out)
	print(file=fh_out)

if args.lineage_file:
	tree = ncbi.get_topology(taxi_id_list)
	fh_linea = open(args.lineage_file, 'w')
	print(tree.get_ascii(attributes=["sci_name", "rank"]), file=fh_linea)
	fh_linea.close()
