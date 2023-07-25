#!/usr/env/python python3
import sys
import matplotlib.pyplot as plt
from icecream import ic
import json

if len(sys.argv) != 2:
    print("Usage: python3 {} data_stat.json".format(sys.argv[0]))
    exit()


def plot_stat(data, title, output):
    plt.figure(figsize=(10,5))#设置画布的尺寸
    plt.title(title,fontsize=20)#标题，并设定字号大小
    plt.boxplot(data)
    plt.savefig(output) #save as png

def plot_hist(data, title, output):
    plt.figure(figsize=(10,5))#设置画布的尺寸
    plt.title(title,fontsize=20)#标题，并设定字号大小
    plt.hist(data, bins=30)
    plt.savefig(output) #save as png

def smart_open(file):
    import gzip
    if file.endswith("gz"):
        return gzip.open(file, 'rt')
    else:
        return open(file)

def load_from_db(handle):
    data = json.load(handle)
    return data

def main():
    total_bases = 0
    total_passes = 0
    # load values from json file
    with open(sys.argv[1], 'r') as fh:
        dataset = load_from_db(fh)
    read_bases = dataset['read_length']
    passes = dataset['passes']
    accuracy = dataset['accuracy']
    plot_hist(read_bases, 'HiFi read length distribution', 'read_length')
    plot_hist(passes, 'pass number distribution', "passes_dist")
    plot_stat(accuracy, 'accuracy value distribution', "accuracy")
    for p,l in zip(passes, read_bases):
        total_passes += int(p)
        total_bases += int(l)

    total_ccs = len(read_bases)
    print(f"total ccs = {total_ccs}")
    print(f"total passes = {total_passes}")
    print(f"total bases = {total_bases}")

if __name__ == '__main__':
    main()


