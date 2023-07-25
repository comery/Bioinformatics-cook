import re
import sys
from igraph import *
from collections import defaultdict
if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} *.path *.gfa")

paths = {}
segment_contigs = {}
node_count = 0

# Get contig paths from contigs.paths
with open(sys.argv[1], 'r') as file:
    name = file.readline()
    path = file.readline()

    while name != "" and path != "":
        while ";" in path:
            path = path[:-2]+","+file.readline()

        start = 'NODE_'
        end = '_length_'
        contig_num = str(int(re.search('%s(.*)%s' % (start, end), name).group(1))-1)

        segments = path.rstrip().split(",")

        if contig_num not in paths:
            node_count += 1
            paths[contig_num] = [segments[0], segments[-1]]
        for segment in segments:
            if segment not in segment_contigs:
                segment_contigs[segment] = set([contig_num])
            else:
                segment_contigs[segment].add(contig_num)

        name = file.readline()
        path = file.readline()
links = []
links_map = defaultdict(set)
# Get contig paths from contigs.paths
with open(sys.argv[2], 'r') as file:
    line = file.readline()

    while line != "":

        # Identify lines with link information
        if "L" in line:
            strings = line.split("\t")
            f1, f2 = strings[1]+strings[2], strings[3]+strings[4]
            links_map[f1].add(f2)
            links_map[f2].add(f1)
            links.append(strings[1]+strings[2]+" "+strings[3]+strings[4])
        line = file.readline()
# Create graph
g = Graph()
# Add vertices
g.add_vertices(node_count)
for i in range(len(g.vs)):
    g.vs[i]["id"]= i
    g.vs[i]["label"]= str(i+1)
for i in range(len(paths)):
    segments = paths[str(i)]

    start = segments[0]
    start_rev = ""
    if start.endswith("+"):
        start_rev = start[:-1]+"-"
    else:
        start_rev = start[:-1]+"+"

    end = segments[1]
    end_rev = ""
    if end.endswith("+"):
        end_rev = end[:-1]+"-"
    else:
        end_rev = end[:-1]+"+"

    new_links = []

    if start in links_map:
        new_links.extend(list(links_map[start]))
    if start_rev in links_map:
        new_links.extend(list(links_map[start_rev]))
    if end in links_map:
        new_links.extend(list(links_map[end]))
    if end_rev in links_map:
        new_links.extend(list(links_map[end_rev]))

    for new_link in new_links:
        if new_link in segment_contigs:
            for contig in segment_contigs[new_link]:
                if i!=int(contig):
                    g.add_edge(i,int(contig))

g.simplify(multiple=True, loops=False, combine_edges=None)
