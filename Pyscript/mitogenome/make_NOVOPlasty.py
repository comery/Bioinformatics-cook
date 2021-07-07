#!/usr/env/python python3
import sys
if len(sys.argv) < 2:
    print("Usage: python3 " + sys.argv[0] + " sample")
    exit()

sample = sys.argv[1]
with open(sample + "/config.txt", 'w') as cf:
    text = """
Project:
-----------------------
"""
    text += "Project name          = {}".format(sample)
    text +="""
Type                  = mito
Genome Range          = 13000-20000
K-mer                 = 39
Max memory            = 30
Extended log          = 0
Save assembled reads  = no
Seed Input            = /hwfssz1/ST_DIVERSITY/P18Z10200N0197_Phylogeny/USER/yangchentao/MT10K/jianglu/Popillia_japonica.mito.fa
Reference sequence    =
Variance detection    = no
Heteroplasmy          =
Chloroplast sequence  =
Dataset 1:
-----------------------
Read Length           = 150
Insert size           = 270
Platform              = illumina
Single/Paired         = PE
Combined reads        =
"""
    text += "Forward reads         = /hwfssz1/ST_DIVERSITY/P18Z10200N0197_Phylogeny/USER/yangchentao/MT10K/jianglu/jianglu_6bettles/insect_T760-S01-01-{}_good_1.fq.gz\n".format(sample)
    text += "Reverse reads         = /hwfssz1/ST_DIVERSITY/P18Z10200N0197_Phylogeny/USER/yangchentao/MT10K/jianglu/jianglu_6bettles/insect_T760-S01-01-{}_good_2.fq.gz\n".format(sample)
    text +="""
Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.8
Insert Range strict   = 1.3
"""
    cf.write(text)

with open(sample + "/novoplasty.sh", 'w') as sh:
    cmd = "perl \
        /hwfssz1/ST_DIVERSITY/PUB/USER/mengguanliang/soft/NOVOPlasty2.7.2/NOVOPlasty2.7.2.pl \
        -c config.txt"
    sh.write(cmd)
