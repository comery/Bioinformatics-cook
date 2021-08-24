#!/usr/bin/env python3
import sys
import os
#from icecream import ic

if len(sys.argv) < 2:
    sys.exit(f"print {sys.argv[0]} accession_number/update")


def findpath(subdir):
    # example, GCF_015227675.2_mRatBN7.2_genomic.fna.gz
    dir = '/hwfssz1/pub/database/ftp.ncbi.nih.gov/genomes/all/'
    tmp = subdir.split("_")
    cate = tmp[0]
    dir1 = tmp[1][0:3]
    dir2 = tmp[1][3:6]
    dir3 = tmp[1][6:9]
    search_dir = os.path.join(dir, cate, dir1, dir2, dir3)
    #ic(search_dir)
    find = 0
    for root, dirs, files in os.walk(search_dir):
        for file in files:
            if subdir in file:
                if 'gnomon_rna.fna.gz' in file:
                    abspath = os.path.join(root, file)
                    #print(abspath)
                    find = 1
    if find:
        return abspath
    else:
        return None


def main():
    if 'GCA' in sys.argv[1] or 'GCF' in sys.argv[1]:
        find = findpath(sys.argv[1])
        if find:
            print(find)
        else:
            print("can not find this one")
    else:
        sys.exit("I can only deal with GCA and GCF")

if __name__ == '__main__':
    main()
