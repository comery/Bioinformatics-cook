#!/usr/bin/env python3
import argparse
#import fileinput
#import numpy as np
import pandas as pd

description = '''
Description

    Transform the coordinate.
    Usage:
	python T1_2.py -i ./temporary_dir/test.makeup.agp -a ./temporary_dir/test2.bed -o new.bed
 
'''

parser = argparse.ArgumentParser(description = description, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", required=True, type=str, help="Input agp file")
parser.add_argument("-a", required=True, type=str, help="Input bed file")
parser.add_argument("-o",required=True, type=str, help="Output bed")
args   = parser.parse_args()

D   = {}
Dor = {}
with open(args.i,'r') as fi:
    for line in fi :        
        l  = line.strip().split()
        D.setdefault(l[0],{})
        D[l[0]][(int(l[1]),int(l[2]))] = (int(l[6]),int(l[7]))
        Dor[l[0]] = l[8]

print(D)

'''
def get_k(dict, value):
    return [k for k, v in dict.items() if v == value]
'''

# assume : pat + mat :in same line and orientation
dfi = pd.read_csv(args.a,sep='\t') #bed file
for i in range(len(dfi)):
    p1,p2 = int(dfi.loc[i][1]),int(dfi.loc[i][2]) #coordinate
    m1,m2 = int(dfi.loc[i][4]),int(dfi.loc[i][5])
    for k in (dfi.loc[i][0],dfi.loc[i][3]) :
        for (a,b) in D[k].keys():
            if Dor[k] == '+' : 
                if a <= p2 <= b :
                    p1,p2 = D[k][(a,b)][0]+(p1-a),D[k][(a,b)][0]+(p2-a)
                    continue
                if a <= m2 <= b :
                    m1,m2 = D[k][(a,b)][0]+(m1-a),D[k][(a,b)][0]+(m2-a)
                    continue
            if Dor[k] == '-' :
                if a <= p2 <= b :
                    p1,p2 = D[k][(a,b)][0]-(p1-a),D[k][(a,b)][0]-(p2-a)
                    continue
                if a <= m2 <= b :
                    m1,m2 = D[k][(a,b)][0]-(m1-a),D[k][(a,b)][0]-(m2-a)

    dfi.iloc[[i,1]] = p1
    dfi.iloc[[i,2]] = p2
    dfi.iloc[[i,4]] = m1
    dfi.iloc[[i,5]] = m2

dfi.to_csv(args.o, sep='\t')


'''
    for j,k in zip(D[pat].keys(),D[mat].keys()): 
        if Dor[pat] == '+' or Dor[mat] == '+':
            if j[0] <= p2 <= j[1] :
                p1,p2 = int(D[pat][j][0]+(p1-j[0])),int(D[pat][j][0]+(p2-j[0]))
            if k[0] <= m2 <= k[1] :
                m1,m2 = int(D[mat][k][0]+(m1-k[0])),int(D[mat][k][0]+(m2-k[0]))
        elif Dor[pat] == '-' or Dor[mat] == '-':
            if j[0] <= p2 <= j[1] :
                p1,p2 = int(D[pat][j][0]-(p1-j[0])),int(D[pat][j][0]-(p2-j[0]))
            if k[0] <= m2 <= k[1] :
                m1,m2 = int(D[mat][k][0]-(m1-k[0])),int(D[mat][k][0]-(m2-k[0]))

    dfi.iloc[[i,1]]   = p1
    dfi.iloc[[i,2]]   = p2
    dfi.iloc[[i,4]]   = m1
    dfi.iloc[[i,5]]   = m2
'''    

