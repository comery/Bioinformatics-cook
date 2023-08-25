import os
import sys
import math
import argparse
from matplotlib import pylab
import pylab as plt
import numpy as np
import pandas as pd
import mappy as mp
import matplotlib.patches as patches

def letter_phred(symbols):
    qual = []
    for i in symbols:
        qual.append(ord(i) - 33)
    return qual


def plot_fastq_qualities(filename, ax=None, limit=10000, output='fastq_quality_plot.png'):
    res=[]
    c=0
    for read in mp.fastx_read(filename, read_comment=False):
        score=letter_phred(read[2])
        res.append(score)
        c+=1
        if c>limit:
            break
    df = pd.DataFrame(res)
    l = len(df.T)+1

    if ax==None:
        f,ax=plt.subplots(figsize=(12,5))
    rect = patches.Rectangle((0,0),l,20,linewidth=0,facecolor='r',alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0,20),l,8,linewidth=0,facecolor='yellow',alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0,28),l,12,linewidth=0,facecolor='g',alpha=.4)
    ax.add_patch(rect)
    df.mean().plot(ax=ax,c='black')
    boxprops = dict(linestyle='-', linewidth=1, color='black')
    df.plot(kind='box', ax=ax, grid=False, showfliers=False,
            color=dict(boxes='black',whiskers='black')  )
    ax.set_xticks(np.arange(0, l, 5))
    ax.set_xticklabels(np.arange(0, l, 5))
    ax.set_xlabel('position(bp)')
    ax.set_xlim((0,l))
    ax.set_ylim((0,40))
    ax.set_title('per base sequence quality')
    plt.savefig(output) #save as png


def fastq_to_dataframe(filename, size=1000):
    """Convert fastq to dataframe.
        size: limit to the first reads of total size
        Returns: dataframe with reads
    """
    i=0
    res=[]
    for read in mp.fastx_read(filename, read_comment=False):
        res.append([read[0], read[1]])
        i += 1
        if i > size:
            break
    df = pd.DataFrame(res, columns=['id','seq'])
    df['length'] = len(read[1])
    return df

def normpdf(x, mean, sd):
    """sample a normal distribution at given point"""

    var = float(sd)**2
    denom = (2*math.pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def plot_fastq_gc_content(filename, ax=None, limit=50000, output='fastq_gc_plot.png'):
    from Bio.SeqUtils import GC
    if ax==None:
        f,ax=plt.subplots(figsize=(12,5))
    df = fastq_to_dataframe(filename, size=limit)
    gc = df.seq.apply(lambda x: GC(x))
    gc.hist(ax=ax,bins=150,color='black',grid=False,histtype='step',lw=2)
    ax.set_xlim((0,100))
    x=np.arange(1,100,.1)
    f = [normpdf(i, gc.mean(), gc.std()) for i in x]
    ax2=ax.twinx()
    ax2.plot(x,f)
    ax2.set_ylim(0,max(f))
    ax.set_title('GC content',size=15)
    plt.savefig(output)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(f"python3 {sys.argv[0]} -g -q *.fq.gz")
    else:
        parser = argparse.ArgumentParser()
        parser.add_argument("-g", "--gc", dest="gc_content", action="store_true",
                            help="plot gc content distribution")
        parser.add_argument("-q", "--qual", dest="quality", action="store_true",
                            help="plot read quality distribution")
        parser.add_argument("fastx", type=str,
                            help="input file, fasta or fastq")
        args = parser.parse_args()
        if not args.gc_content and not args.quality:
            sys.exit("what do you want to draw?")
        if args.gc_content:
            plot_fastq_gc_content(args.fastx, ax=None, limit=50000)
        if args.quality:
            plot_fastq_qualities(args.fastx, ax=None, limit=50000)
