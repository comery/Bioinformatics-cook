#!/usr/bin/env python

from posixpath import realpath
import sys
import os
import argparse
from ete3 import Tree, PhyloTree, TreeStyle, TextFace, NodeStyle # install ete3 with PyQt5
from Bio import SeqIO
from Bio.Seq import Seq
from icecream import ic


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

def add2dim(thedict, ka, kb, val):
    if ka in thedict:
        thedict[ka].update({kb: val})
    else:
        thedict.update({ka: {kb: val}})


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
    start_str = "Bayes Empirical Bayes (BEB) analysis"
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
        line = line.strip()
        if line:
            codon_n, aa, p = line.rstrip("\n").split(" ")
            codon_n = int(codon_n)
            p = float(p.rstrip("*"))
            #if aa != '-':
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


def trans_aln(aln_region):
    # input an coding aln
    aa_aln_str = ""
    for seq_obj in aln_region:
        dna = Seq(str(seq_obj.seq))
        aa_s = dna.translate(gap="-")
        aa_aln_str += '>{}\n{}\n'.format(seq_obj.id, aa_s)
    return aa_aln_str

def convertPos(aln, ps_sites):
    # get all positive sites
    pos = []
    for item in ps_sites:
        pos.append(item[0])
    RealPosMap = {}
    Mutations = {}
    for record in aln:
        species = record.id
        seq = record.seq
        protein = seq.translate(gap="-")   # nucleotide to amino acid
        for p in pos:
            tmp = protein[:p-1]
            indel = tmp.count('-')
            if indel == len(tmp) and protein[p-1] == "-":  # the sequence is like '---------' at starting
                real_pos = -1
            else:
                real_pos = p - indel
            amino_acid = protein[p-1]
            add2dim(RealPosMap, p, species, real_pos)
            add2dim(Mutations, p, species, amino_acid)
    return RealPosMap, Mutations

def say_real_pos(codon, RealPosMap, Mutations):
    out = ""
    for sp in RealPosMap[codon]:
        out += f"{sp}:{RealPosMap[codon][sp]}:{Mutations[codon][sp]}" + ";"
    return out


def draw_tre_and_aln(aln, tre_str, ps_site, og, target_sp, extend_aa):
    """
    draw tree and alignment together
    """
    site, aa, score=ps_site

    # get leaves from tree string
    tstax = getleaf(tre_str)

    ts = TreeStyle()
    ts.margin_left = 5
    ts.margin_right = 30
    ts.margin_top = 20
    ts.tree_width = 50

    codon_start = site * 3 - 2
    codon_end = site * 3
    # up and down extend extend_aa aa 0 based or 1 based
    region_start = codon_start - extend_aa*3 - 1
    region_end = codon_end + extend_aa*3
    # fix if region_start exceeds the start issue
    if region_start < 0:
        region_start = 0
    aln_region = aln[:,region_start:region_end]
    #print("Showing the [codon: {}] of alignment region from {} to {}".format(site, region_start+1, region_end))
    #aln_region_str = aln_region.format("fasta") # remove this useage in biopython latest version
    aln_region_str = format(aln_region, "fasta")
    t = PhyloTree(tre_str, alignment=aln_region_str, alg_format="fasta")

    # interfact
    def _set_style(t):
        # input an t with alignment, add label to it
        info = TextFace("{}\nCodon:{}\nScore:{}".format(og,site,score), fsize=8, fgcolor='black', ftype='Arial')
        info.margin_top = 10
        info.margin_right = 20
        info.margin_left = 5
        t.add_face(info, column=0, position="branch-bottom")
        #t.add_face(TextFace("Codon:{}".format(site)),column=0,position="branch-bottom")
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
    aa_aln_region_str = trans_aln(aln_region)
    t_aa = PhyloTree(tre_str, alignment=aa_aln_region_str, alg_format="fasta")
    t_aa = _set_style(t_aa)

    return t, t_aa, ts

def check_alignment(aln, codon):
    aln_aa = trans_aln(aln).strip()
    aas = []
    for record in aln_aa.split("\n"):
        if record.startswith(">"):
            continue
        else:
            aas.append(record[codon-1])
    gap_count = aas.count('-')
    aas_set = set(aas)
    aa_categories = len(list(aas_set))
    return gap_count, aa_categories

def main():
    mlc_f = args.input
    aln_f = args.fasta
    subdir = args.prefix
    target_sp = args.foreground
    outdir = os.path.join(args.outdir, subdir)
    outdir = os.path.abspath(outdir)

    if os.path.exists(outdir) == False:
        os.makedirs(outdir)

    for d in ['filter', 'pass']:
        tmp_d = os.path.join(outdir, d)
        if os.path.exists(tmp_d) == False:
            os.makedirs(tmp_d)

    # write tree to file
    tree_file = os.path.join(outdir, subdir + ".tree")
    aln_obj, tre_str = get_aln_tre_sites(mlc_f, aln_f)
    record_tree(tre_str, tree_file)

    # write positive sites to file
    positive_sites_file = os.path.join(outdir, subdir + ".sites")
    ps = open(positive_sites_file, 'w')
    psg_sites = parse_BEB_sites(mlc_f)

    # get realposmap
    RealPosMap, Mutations = convertPos(aln_obj, psg_sites)

    if len(psg_sites) == 0:
        print("{} has no sites".format(subdir))
        exit(1)

    for psg_site in psg_sites:
        codon, aa, score = psg_site
        ic(codon, aa, score)
        gap_count, aa_categories = check_alignment(aln_obj, codon)
        position = say_real_pos(codon, RealPosMap, Mutations)
        # output site quality
        tmp_outdir = ''
        if gap_count > args.gap or aa_categories > args.ka:
            ps.write(f"{subdir}\t{codon}\t{position}\t{aa}\t{score}\tfilter\t{gap_count}\t{aa_categories}\n")
            tmp_outdir = 'filter'
        else:
            ps.write(f"{subdir}\t{codon}\t{position}\t{aa}\t{score}\tpass\t{gap_count}\t{aa_categories}\n")
            tmp_outdir = 'pass'

        out_png = "{famid}_Codon{c}_{p}.pdf".format(famid=subdir,c=codon,p=score)
        out_aa_png = "{famid}_Codon{c}_{p}_AA.pdf".format(famid=subdir,c=codon,p=score)
        out_png_pth = os.path.join(outdir, tmp_outdir, out_png)
        out_aa_png_pth = os.path.join(outdir, tmp_outdir, out_aa_png)
        t, t_aa, ts = draw_tre_and_aln(aln_obj, tre_str,
                                        psg_site, subdir,
                                        target_sp, args.extend)
        t.render(out_png_pth, w=args.width, units="mm", dpi=args.dpi)
        t_aa.render(out_aa_png_pth, w=args.width, units="mm", dpi=args.dpi)

    ps.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='input an og draw the positively selected sites for manual check of reliable sites and alignment issue sites')
    parser.add_argument('input', type=str, help="H1.mlc file generated from PAML")
    parser.add_argument('fasta', type=str, help="cds sequence alignment with fasta format")
    parser.add_argument('foreground', type=str, help="specific short name for #1")
    parser.add_argument('prefix', type=str, help='subdir for ortholog')
    parser.add_argument('-o', '--outdir', type=str, dest='outdir', default="Results_phyaln", required=False, help='outdir, default=Results_phyaln')
    parser.add_argument('-e', '--extend', type=int, dest='extend', default=8, help='extended amino acid length for image, default=8')
    parser.add_argument('-g', '--gap', type=int, dest='gap', default=1, help="max gap number allowed in a site, default=1")
    parser.add_argument('-k', '--ka', type=int, dest='ka', default=2, help="max amino acid catelogies allowed in a site, default=2")
    parser.add_argument('-w', '--width', type=int, dest='width', default=180, help="image width, default=180mm")
    parser.add_argument('-d', '--dpi', type=int, dest='dpi', default=300, help="image dpi, default=300")


    args = parser.parse_args()
    main()
