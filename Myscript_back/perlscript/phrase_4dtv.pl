#!/usr/bin/perl -w
# Copyright (c) njau 2011
# Writer:         songxm

my $ver="1.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作
my %opts;
GetOptions(\%opts,"i1=s","o=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i1}) ||!defined($opts{o}) || defined($opts{h}))
{
 print <<" Usage End.";
 Description:
  
  Version: $ver

 Usage:

  -i1  infile1    must be given
  -o  outfile    must be given
  
  -h   Help document

 Usage End.

 exit;
}
############# perl four_ssite.pl -i1 dna.fa -o dan.fa.out
#####提取DNA序列中符合四倍简并位点的序列；
my $programe_dir=basename($0);
my $infile1=$opts{i1};
my $outfile=$opts{o};
my @bases;my %hash1;my %hash2;my %hash3;my %hash4;

open (IN1,"$infile1") || die "Can't open $infile1\n";
open (OUT2,">$outfile") ||die "can't open the $outfile\n";
$/=">";
<IN1>;
while (<IN1>)
{
 chomp;
 my ($name,$seq)=(split(/\n/,$_))[0,1];
 my %hash=(
  'TC'=>'Ser',
  'CT'=>'Leu',
  'CC'=>'Pho',
  'CG'=>'Arg',
  'AC'=>'Thr',
  'GT'=>'Val',
  'GC'=>'Phe',
  'GG'=>'Gln',
  );
 my $len=length $seq;
 for (my $i=0;$i<$len ;$i+=3)
 {
  $seq=~m/(\w\w)(\w)/g;

 if (defined $hash{$1})
 {
  print "$1\t$2\n";
  $hash2{$name}.=$2;
 }
 }

}

close IN1;

foreach  (keys %hash2)
{
 print OUT2 ">$_\n$hash2{$_}\n";

}
close OUT2;