from Bio import Phylo
from io import StringIO
tree = Phylo.read('tree.nwk', 'newick')
#Phylo.draw_ascii(tree)
#for leaf in tree.get_terminals():
 #   print(leaf)

tre_str = "(COW_R00552_COW: 0.303352, (MOUSE_R07609_MOUSE: 1.041145, ((((MARMOS_R17819_MARMOS: 0.595782, AOTNAN_R19963_AOTNAN: 0.050694): 0.002406, (SAIBOL_R05912_SAIBOL: 0.063089,CEBIMI_R03570_CEBIMI: 0.045699): 0.008387): 0.069147, (MACACA_R04093_MACACA: 0.037942, (CHIMPA_R19344_CHIMPA: 0.009620, HUMAN_R12548_HUMAN: 0.016771): 0.034695): 0.034348): 0.211448,TRESHR_R00513_TRESHR: 0.292816): 0.093355): 0.279026);"
tree1 = Phylo.read(StringIO(tre_str), "newick")
#for l in tree1.get_terminals():
#    print(i)

def getleaf(tre_str):
    tmp = tre_str.split(",")
    for i in tmp:
        #print(i)
        i = i.replace("(", "").replace(")", "").replace(":", "")
        if  "_" in i:
            if " " in i:
                print(i.strip().split()[0])
            else:
                print(i)

getleaf(tre_str)
