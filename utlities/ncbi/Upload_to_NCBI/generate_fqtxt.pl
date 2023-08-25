#!/usr/bin/perl
use strict;

die "perl generate_fqtxt.pl <library> <fqlist>\n" if @ARGV<2;

my $library=shift;
my $fqlist=shift;
my %lib;

open LI,$library;
while(<LI>){
	chomp;
	my @A=split /\s+/;
	$lib{$A[1]}=[$A[0],$A[2]];
}
close LI;

my $head="#sampleID\tBatch\tDate\tLib\tInsertSize\t-SD/+SD\tlength(read1/read2)\tlane\tevaluation\tfqSize\tfqPix\tRead_Num\tFq1_md5\tFq2_md5\n";
print $head;

open FQ,$fqlist;
while(<FQ>){
	chomp;
	my $file=(split /\//)[-1];
	next unless $file=~/(\d+)_(I\d+)_(FC\w+)_(L\d+)_(.+)_1\.fq\.gz$/;
	my $date=$1; my $libid=$5; my $smpid=$lib{$libid}->[0]; my $libsize=$lib{$libid}->[1];
	my $sdsize;
	$sdsize='-10/+10' if $libsize < 170;
	$sdsize='-10/+10' if $libsize==170;
	$sdsize='-20/+20' if $libsize==200;
	$sdsize='-20/+10' if $libsize==250;
	$sdsize='-20/+20' if $libsize==500;
	$sdsize='-50/+50' if $libsize==800;
	$sdsize='-200/+200' if $libsize==2000;
	$sdsize='-500/+500' if $libsize==5000;
	$sdsize='-2000/+2000' if $libsize==10000;
	$sdsize='-5000/+5000' if $libsize==20000;
	$sdsize='-10000/+10000' if $libsize==40000;
	my $prefix=$1 if /^(.+)_1\.fq\.gz$/;
	my $lane=(split /\//,$prefix)[-1];
	my $md51=`md5sum $_|cut -d " " -f1`;
	chomp $md51;
	my $fqsize1=`ls -lh $_|cut -d " " -f5`;
	chomp $fqsize1;
	my $fq2=$prefix.'_2.fq.gz';
	my $md52=`md5sum $fq2|cut -d " " -f1`;
	chomp $md52;
	my $fqsize2=`ls -lh $fq2|cut -d " " -f5`;
	chomp $fqsize2;
	my $fqsize=$fqsize1.'+'.$fqsize2;
	my $readlen = `gzip -dc $fq2|head -2 |tail -1 |wc -m `;
	$readlen=~s/\s+//g;
	my $readnum = (split /\s+/,`gzip -dc $fq2|wc -l  `)[0];	
	$readnum/=2;
	$readnum=$readnum*$readlen;
	$readlen="$readlen/$readlen";
	print "$smpid\t1\t$date\t$libid\t$libsize\t$sdsize\t$readlen\t$lane\t2NE\t$fqsize\t$prefix\t$readnum\t$md51\t$md52\n";
}
close FQ;

