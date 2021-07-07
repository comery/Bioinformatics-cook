#!/usr/env/python python3
import argparse
from Bio import Phylo
import palettable.cartocolors.qualitative

try:
    from palettable.tableau import GreenOrange_12
except ImportError:
    print("We need palettable, sorry...\ninstall: pip install palettable")
    sys.exit(1)

try:
    import palettable.cartocolors.qualitative
except ImportError:
    print("We need palettable, sorry...\ninstall: pip install palettable")
    sys.exit(1)

parser = argparse.ArgumentParser(
    description="automatically generate a color configure for itol based on\n"
    + "your tree leaves, for example, leaf name is supposed to be XX_XX_family",
)

parser.add_argument(
    "-tree",
    metavar="FILE",
    type=str,
    required=True,
    help="phylogenetic tree file.",
)
parser.add_argument(
    "-sep",
    metavar="STR",
    type=str,
    required=False,
    default="_",
    help="separator in leaf name",
)
parser.add_argument(
    "-index",
    metavar="INT",
    type=int,
    required=False,
    default=3,
    help="general configures to circos parameters",
)
parser.add_argument(
    "-type",
    metavar="STR",
    type=str,
    required=True,
    help="which type configure do you want to generate?"
    + "[range|label|label_background]",
)


args = parser.parse_args()


""" colors theme """
colors = []
print(
    "colors theme is designed by palettable." +\
    " \nsee more:http://colorbrewer2.org\n")

colors = palettable.cartocolors.qualitative.Antique_10.colors
colors += palettable.cartocolors.qualitative.Bold_10.colors
colors += palettable.cartocolors.qualitative.Pastel_10.colors
colors += palettable.cartocolors.qualitative.Prism_10.colors
colors += palettable.cartocolors.qualitative.Safe_10.colors
colors += palettable.cartocolors.qualitative.Vivid_10.colors

def rgb2hex(rgb):
    r, g, b = int(rgb[0]), int(rgb[1]), int(rgb[2])
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)



leaves = []
families = []
tree = Phylo.read(args.tree, 'newick')
#Phylo.draw_ascii(tree)
for leaf in tree.get_terminals():
    tmp = str(leaf).split(args.sep)
    leaves.append(str(leaf))
    if len(families) == 0:
        families.append(tmp[args.index - 1])
    elif tmp[args.index - 1] not in families:
        families.append(tmp[args.index - 1])

if len(families) > len(colors):
    print("Too many families to generate enough colors")
    exit()

# map colors to families 
design = {}
for i in range(len(families)):
    design[families[i]] = colors[i]

# it must be .txt file output
with open(args.tree + ".color.conf.txt", 'w') as fh:
    fh.write("TREE_COLORS\nSEPARATOR SPACE\nDATA\n\n")
    for e in leaves:
        tmp = e.split(args.sep)
        co = design[tmp[args.index - 1]]
        #this_color = ",".join(str(v) for v in co)
        this_color = rgb2hex(co)
        fh.write("{0} {1} {2}\n".format(e, args.type, this_color))

print("All done !")
