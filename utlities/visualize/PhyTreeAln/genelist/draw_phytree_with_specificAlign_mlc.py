#!/usr/bin/env python

# input an og draw the positively selected sites for manual check of reliable sites and alignment issue sites
import sys
import os
from ete3 import Tree, PhyloTree, TreeStyle, TextFace, NodeStyle # install ete3 with PyQt5
from Bio import SeqIO
from Bio.Seq import Seq

# parse the mlc and return the 
# alignment and the slection tre and sites

def get_seq_mlc(og_path):
    family_id = og_path.split("/")[-1]
    mlc_f = os.path.join(og_path, "H1.mlc")
    aln_f = os.path.join(og_path, family_id + ".seq") # this is phylip fromat
    return family_id, mlc_f, aln_f


def getleaf(tre_str):
    """
    get leaves from tree string
    actually you can use
    Phylo.read(StringIO(tre_str), "newick")
    but it seems something wrong!
    """
    leaves = []
    tmp = tre_str.split(",")
    for i in tmp:
        i = i.replace("(", "").replace(")", "").replace(":", "")
        if  "_" in i:
            if " " in i:
                leaves.append(i.strip().split()[0])
            else:
                leaves.append(i)
    return leaves

def record_tree(string, outfile):
    with open(outfile, 'w') as fh:
        fh.write(string)

def phylip2fasta(inpath, outfile):
    try:
        infile = open(inpath,"r")
    except OSError:
        exit("There was a problem opening the specified directory. Are you sure it exists? Quitting.")

    out = open(outfile,"w")

    firstline = infile.readline()
    while True:
        name = infile.readline().strip()
        seq = infile.readline().strip()
        if name:
            out.write(">" + name + "\n")
            out.write(seq + "\n")
        else:
            break

    infile.close()
    out.close()

def parse_BEB_sites(mlc_f):
    """
    parse positive sites from mlc file
    """
    with open(mlc_f, 'r') as fh:
        lines = fh.readlines()
    # filter lines between 
    # [Positive sites for foreground lineages Prob(w>1):]
    # and 
    # [The grid (see ternary graph for p0-p1)]
    # start_str = "Positive sites for foreground lineages Prob(w>1):"
    psg_sites = []
    start_str = "Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)"
    end_str = "The grid (see ternary graph for p0-p1)"
    start_line = 0
    end_line = 0
    for i, line in enumerate(lines):
        if line.startswith(start_str):
            # also NEB will match
            # use previous line
            start_line = i+2
        elif line.startswith(end_str):
            end_line = i-1
    for line in lines[start_line:end_line]:
        line = line.lstrip().rstrip()
        if line:
            codon_n, aa, p = line.rstrip("\n").split(" ")
            codon_n = int(codon_n)
            p = float(p.rstrip("*"))
            if aa != '-':
                psg_sites.append([codon_n, aa, p])

    return sorted(psg_sites)


def get_aln_tre_sites(mlc_f,aln_f):
    """
    input an mlc file and it's alnignment file
    return the alignment (biopython object) and tree string
    # the '-' in the first sequence were omited
    """
    import os
    from Bio import AlignIO
    import re
    from Bio.Phylo.PAML import codeml

    og = mlc_f.split("/")[-2]
    aln_obj = AlignIO.read(aln_f, "fasta")

    ## Tree
    mlc_text = ''
    with open(mlc_f, 'r') as f:
        mlc_text = f.read()
    mlc = codeml.read(mlc_f)
    try:
        tre_str = mlc['NSsites'][2]['tree']
    except KeyError as e:
        tre_str = ''
        print("No phylotree in the codeml mlc out file!")

    return aln_obj, tre_str

def draw_tre_and_aln(aln_region, tre_str, og, target_sp):
    """
    draw tree and alignment together
    """
    #site, aa, score=ps_site

    # get leaves from tree string
    tstax = getleaf(tre_str)

    ts = TreeStyle()
    ts.margin_left = 5
    ts.margin_right = 30
    ts.margin_top = 20
    ts.tree_width = 50
    """
    codon_start = site * 3 - 2
    codon_end = site * 3
    # up and down extend 5aa 0 based or 1 based
    region_start = codon_start - 5*3 - 1
    region_end = codon_end + 5*3
    # fix if region_start exceeds the start issue
    if region_start < 0:
        region_start = 0
    aln_region = aln[:,region_start:region_end]
    """
    #print("Showing the [codon: {}] of alignment region from {} to {}".format(site, region_start+1, region_end))
    aln_region_str = aln_region.format("fasta")

    t = PhyloTree(tre_str, alignment=aln_region_str, alg_format="fasta")

    # interfact
    def _set_style(t):
        # input an t with alignment, add label to it
        """
        info = TextFace("{}\nCodon:{}\nScore:{}".format(og,site,score), fsize=8, fgcolor='black', ftype='Arial')
        info.margin_top = 10
        info.margin_right = 20
        info.margin_left = 5
        t.add_face(info, column=0, position="branch-bottom")
        #t.add_face(TextFace("Codon:{}".format(site)),column=0,position="branch-bottom")
        """
        ## label the longbranch
        nstyle = NodeStyle()
        # red line
        #nstyle["bgcolor"] = "DarkSeaGreen"
        #nstyle["bgcolor"] = "LightSalmon"
        nstyle["hz_line_type"] = 0
        #nstyle["hz_line_color"] = "#ff0000"
        for tst in tstax:
            tsnode = t.get_leaves_by_name(name=tst)[0]
            tsnode.set_style(nstyle)
            # add #1 to target species
            if target_sp in tsnode.name:
                tsnode.name = tsnode.name + "_#1"
        return t
    t = _set_style(t)
    ## add AA alignment
    def _trans_aln(aln_region):
        # input an coding aln
        aa_aln_str = ""
        for seq_obj in aln_region:
            dna = Seq(str(seq_obj.seq))
            aa_s = dna.translate(gap="-")
            aa_aln_str += '>{}\n{}\n'.format(seq_obj.id, aa_s)
        #print(aa_aln_str)
        return aa_aln_str

    aa_aln_region_str = _trans_aln(aln_region)
    t_aa = PhyloTree(tre_str, alignment=aa_aln_region_str, alg_format="fasta")
    t_aa = _set_style(t_aa)

    return t, t_aa, ts


def main():
    with open(ogs_list, 'r') as fh:
        for i in fh:
            og = i.strip()
            og_path = os.path.abspath(og)
            family_id, mlc_f, aln_p = get_seq_mlc(og_path)
            tmp_outdir = os.path.join(outdir, family_id)
            if os.path.exists(tmp_outdir) == False:
                os.makedirs(tmp_outdir)
            aln_f = os.path.join(tmp_outdir, family_id + ".fasta")
            tree_file = os.path.join(tmp_outdir, family_id + ".tree")
            phylip2fasta(aln_p, aln_f) # convert phylip to fasta
            aln_obj, tre_str = get_aln_tre_sites(mlc_f, aln_f)
            record_tree(tre_str, tree_file) # write tree to file
            """
            for psg_site in psg_sites:
                codon, aa, score = psg_site
                ps.write("{}\t{}\t{}\t{}\n".format(family_id, codon, aa, score))
                out_png = "{famid}_Codon{c}_{p}.png".format(famid=family_id,c=codon,p=score)
                out_aa_png = "{famid}_Codon{c}_{p}_AA.png".format(famid=family_id,c=codon,p=score)
                out_png_pth = os.path.join(tmp_outdir, out_png)
                out_aa_png_pth = os.path.join(tmp_outdir, out_aa_png)
                t, t_aa, ts = draw_tre_and_aln(aln_obj, tre_str,
                                               psg_site, family_id,
                                               target_sp)
                t.render(out_png_pth, w=183, units="mm", dpi=300)
                t_aa.render(out_aa_png_pth, w=183, units="mm", dpi=300)
            """
            codon_start = site * 3 - 2
            codon_end = site * 3
            # up and down extend 5aa 0 based or 1 based
            region_start = codon_start - draw_lenth*3 - 1
            region_end = codon_end + draw_lenth*3
            # fix if region_start exceeds the start issue
            if region_start < 0:
                region_start = 0
            aln_region = aln_obj[:,region_start:region_end]
            out_png = "{famid}_{site}_{start}_{end}.png".format(famid=family_id,
                                                                site=site,
                                                                start=region_start,
                                                                end=region_end)
            out_aa_png = "{famid}_{site}_{start}_{end}_AA.png".format(famid=family_id,
                                                                      site=site,
                                                                      start=region_start,
                                                                      end=region_end)
            out_png_pth = os.path.join(tmp_outdir, out_png)
            out_aa_png_pth = os.path.join(tmp_outdir, out_aa_png)
            t, t_aa, ts = draw_tre_and_aln(aln_region, tre_str,
                                           family_id, target_sp)
            t.render(out_png_pth, w=183, units="mm", dpi=300)
            t_aa.render(out_aa_png_pth, w=183, units="mm", dpi=300)


            ps.close()

if __name__ == "__main__":
    usage = """
    python3 {} <gene_list> <specific short name for #1> <draw_pos> <radius> <outdir|default=Results_phyaln>
    - draw_pos is targe location of codon (aa) in alignment, not nucleotide
    - radius defines the range length.
    """.format(sys.argv[0])
    if len(sys.argv) < 5:
        sys.exit(usage)
    elif len(sys.argv) > 5:
        outdir = sys.argv[5]
        if os.path.exists(outdir) == True:
            print("{} has already existed, be careful!".format(outdir))
    else:
        outdir = "Results_phyaln"
    ogs_list = sys.argv[1]
    target_sp = sys.argv[2]
    site = int(sys.argv[3])
    draw_lenth = int(sys.argv[4])
    main()
