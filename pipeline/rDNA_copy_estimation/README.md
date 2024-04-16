### Estimation rDNA copy number from NGS data using kmer-based method

```shell
# generate all kmers in NGS data
date
export PATH=/your_path/software/Genome/KMC/bin:$PATH
[ -d tmp  ] || mkdir tmp
echo "/your_path/00.dataset/illumina/clean/MF2/D2120626A_L4_110A10.R1.clean.fq.gz
/your_path/00.dataset/illumina/clean/MF2/D2120626A_L4_110A10.R2.clean.fq.gz" > FILES

#计算31 mer k-mer频率
kmc -k31 -t16 -m12 -ci1 -cs100000 @FILES kmcdb tmp 
#生成k-mer频数直方表sample.histo和k-mer直方图
kmc_tools transform kmcdb histogram sample.histo -cx100000 
kmc_dump -ci2 -cx100000  kmcdb ngs.31mer.freq.txt

# 将rDNA ref序列打断成 31mer
python3 /your_path/bin/generateKmer.py ../human_rdna.ref.fa 31 > rdnaRef.31mer.txt
python3 /your_path/bin/extractKmerFreqFromDb.py ngs.31mer.freq.txt rdnaRef.31mer.txt >rdnaRef.31mer.freq_in_rawreads.txt

# for plot
awk '{if ($3/64>500) print $1"\t500"; else print $1"\t"$3/44}' rdnaRef.31mer.freq_in_rawreads.txt > plot.txt
python3 /your_path/bin/visualize/estimate_rDNA_copy.py  rdnaRef.31mer.freq_in_rawreads.txt rdnaRef.31mer.copy_est.pdf 44
date

# control gene, TBP1， 这是一个human的单拷贝基因作为control
python3 /your_path/bin/generateKmer.py ../TBP1_human_seq.dnas 31 > TBP1.31mer.txt
python3 /your_path/bin/extractKmerFreqFromDb.py ngs.31mer.freq.txt >TBP1.31mer.freq_in_rawreads.txt

# plot
awk '{if ($3/64>500) print $1"\t500"; else print $1"\t"$3/44}' TBP1.31mer.freq_in_rawreads.txt > TBP1.plot.txt
python3 /your_path/bin/visualize/estimate_rDNA_copy.py  TBP1.31mer.freq_in_rawreads.txt TBP1.31mer.copy_est.pdf 44

```

