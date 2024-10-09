
[ -d after ] || mkdir after
[ -d before ] || mkdir before

cd before
reads=xxx
ref=xxx
# https://github.com/PacificBiosciences/hg002-ccs/tree/master/concordance
bin=xxx
num_threads=8
minimap2 -a -Q --eqx --secondary=no  -K4G -t $num_threads \
	$ref $reads | samtools view -Sb \
	-@$num_threads > mappings.bam

samtools sort -@$num_threads mappings.bam -o mappings.sorted.bam
samtools index -@$num_threads mappings.sorted.bam

$bin/bamConcordance $ref mappings.sorted.bam out.csv
awk '{print $0",before"}' out.csv > out.before.csv

cd ../after

cd before
reads=xxx
ref=xxx
num_threads=8
minimap2 -a -Q --eqx --secondary=no  -K4G -t $num_threads \
	$ref $reads | samtools view -Sb \
	-@$num_threads > mappings.bam

samtools sort -@$num_threads mappings.bam -o mappings.sorted.bam
samtools index -@$num_threads mappings.sorted.bam

$bin/bamConcordance $ref mappings.sorted.bam out.csv
awk '{print $0",after"}' out.csv > out.after.csv

cd ..

head -1 before/out.before.csv|sed 's/before/is_corr/g' > header.txt
cat before/out.before.csv after/out.after.csv |grep -v '^#' >tmp
cat header.txt tmp > out.csv

Rscript plot.R

