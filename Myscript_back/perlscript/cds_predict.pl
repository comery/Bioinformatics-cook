#!/usr/bin/perl -w
use strict;
use File::Basename;
die "perl $0 <*.fa> <dir of xls> <outdir>" unless (@ARGV == 3);

my $fasta = $ARGV[0];
my $xls = $ARGV[1];
my $outdir = $ARGV[2];
my $fa = basename($fasta);
my $ESTScan = "/ifs4/NGB_ENV/USER/yangchentao/software/RNA-tools/ESTScan";


open OUT, ">cds_predict.sh" ;

print OUT <<CDS
echo ==========start at : `date` ==========
perl /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/CDSpredict/get_cds_blast.pl -fa $fasta -xls $xls/All-Unigene.fa.blast.nr.xls,$xls/All-Unigene.fa.blast.swissprot.xls,$xls/All-Unigene.fa.blast.kegg.xls,$xls/All-Unigene.fa.blast.cog.xls -out $outdir/Blast/All-Unigene.blast -L 100
if [ -d $outdir/ESTscan/estscan ];then rm -rf $outdir/ESTscan/estscan;fi
mkdir -p $outdir/ESTscan/estscan
mv $outdir/Blast/All-Unigene.blast.mrna.fa $outdir/ESTscan/estscan/mrna.seq
perl $ESTScan/prepare_data -e $outdir/ESTscan.conf
perl $ESTScan/build_model $outdir/ESTscan.conf
$ESTScan/estscan $outdir/Blast/All-Unigene.blast.no.fa -o $outdir/ESTscan/All-Unigene.ESTscan.cds.fa.score -t $outdir/ESTscan/All-Unigene.ESTscan.protein.fa.score -M $outdir/ESTscan/estscan/Matrices/*.smat
perl $ESTScan/clear.score.pl $outdir/ESTscan/All-Unigene.ESTscan.cds.fa.score $outdir/ESTscan/All-Unigene.ESTscan.cds.fa -debug
perl $ESTScan/clear.score.pl $outdir/ESTscan/All-Unigene.ESTscan.protein.fa.score $outdir/ESTscan/All-Unigene.ESTscan.protein.fa -debug
cat $outdir/Blast/All-Unigene.blast.cds.fa ESTscan/All-Unigene.ESTscan.cds.fa >$outdir/all.cds.fa
perl /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/CDSpredict/../Denovo/fa_quality.pl -len -Head -N -gc $outdir/all.cds.fa
perl /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/CDSpredict/../Denovo/barplot.pl $outdir/all.cds.fa.quality.xls All-Unigene.cds
perl /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/CDSpredict/cds_stat.pl -blast $outdir/Blast/All-Unigene.blast.cds.fa -estscan $outdir/ESTscan/All-Unigene.ESTscan.cds.fa -output $outdir/PredictSummary.xls

echo ==========start at : `date` ==========
echo Still_waters_run_deep >cds_predict.sh.sign
CDS

close OUT

