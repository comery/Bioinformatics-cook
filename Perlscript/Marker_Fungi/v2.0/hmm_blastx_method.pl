#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use File::Basename;
use Cwd 'abs_path';
die "perl $0 <target.fa> <hmm_profile.lst> <reference.fa[protein]> <outdir>" unless (@ARGV==4);
my $fasta = shift;
my $hmmf = shift;
my $ref = shift;
my $outdir = shift;
$fasta = abs_path($fasta);
$hmmf = abs_path($hmmf);
$ref = abs_path($ref);
$outdir = abs_path($outdir);
my $fa = basename($fasta);
my $qpara = "--queue st.q -P P18Z10200N0197";
#------------------------------------------

`mkdir $outdir` unless ( -d "$outdir");

open OUT, ">hmm_blastx_method.sh";

print OUT <<HMM;
echo ==========start at : `date` ==========
#step 1  hmmserach
#perl $Bin/rename.pl $fasta >$outdir/$fa
[ -d $outdir/hmm_tmp ] && rm -rf $outdir/hmm_tmp
mkdir -p $outdir/hmm_tmp
[ -f hmm.sh ] && rm hmm.sh
cat $hmmf|while read a;do b=`basename \$a`;echo "hmmsearch --cpu 2 --domtblout $outdir/hmm_tmp/\$b.res \$a $outdir/$fa " >>hmm.sh ;done
#sh hmm.sh
perl $Bin/../qsub-sge.pl $qpara --resource vf=0.5G,p=5 -lines 400 -maxjob 50 hmm.sh
find $outdir/hmm_tmp -name "*.res" |while read a;do grep -v "#" \$a |head -1|awk '{print \$1,\$4,\$18,\$19}'|sed "s/.out//" ;done >$outdir/hmm.all.results
awk '{print \$1}' $outdir/hmm.all.results >$outdir/hmm.best.id
perl $Bin/fishInWinter.pl -bf table -ff fasta -gene $outdir/hmm.best.id $outdir/$fa >$outdir/hmm.best.fa
echo "hmmsearch done!"
#step 2 blastx
[ -f "$ref.pin" ] || formatdb -i $ref -o F -p T
echo "blastall -p blastx -e 1e-5 -i $outdir/hmm.best.fa -d $ref -F F -a 4 -o $outdir/blastx.out -m 8" >blastx_vs_ref.sh
perl $Bin/../qsub-sge.pl $qpara --resource vf=1G,p=4 blastx_vs_ref.sh
echo "blastx done!"
# step 3 compare hmm result and blastx's
perl $Bin/best_from_blast_out.pl $outdir/blastx.out >$outdir/blastx.best.out
perl $Bin/best_to_best.pl $outdir/hmm.all.results $outdir/blastx.best.out >$outdir/best_to_best.id
echo "compare done!"
echo ==========start at : `date` ==========
HMM

close OUT;


