#!/usr/bin/perl
#!/bin/bash
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
#use lib "/home/chuanjun/usr/lib/perl5/site_perl/5.8.5";
use Compress::Zlib;

sub usage {
	print <<USAGE;
usage:
	perl $0 [options]
author
	nixiaoming   nixiaoming\@genomics.cn
description:

options:
	-help: print help info
	-fq1    *   <str>  : input 1.fq
	-fq2    *   <str>  : input 2.fq
	-out    *   <str>  : output prefix
	                     out_1.fq out_2.fq 
	                     out.dup.list   out.dup.stat 
	-Rate   *   <str>  : the minimum rate between the num of after_dup_reads && befor_dup_reads ####added by dushuo
e.g.:
	perl $0 -fq1 input_1.fq -fq2 input_2.fq -out  output
USAGE
}

#my ($help,$fq1,$fq2,$out);
my ($help,$fq1,$fq2,$out,$Rate);###added by dushuo
GetOptions(
	"help"=>\$help,
	"fq1=s"=>\$fq1,
	"fq2=s"=>\$fq2,
	"out=s"=>\$out,
	"Rate=f"=>\$Rate,       #############added by dushuo
);
if (defined $help || (!defined $fq1) || (!defined $fq2) || (!defined $out)) {
	&usage();
	exit 0;
}

system "mkdir -p $out && rmdir $out";
open OUT,'>',"$out.dup.stat" or die "can't open $out.dup.stat \nDied ";
print OUT "Start time:".time."\n";

$Rate   ||=0.7;     #######added by dushuo######can be changed based on the quality of PCR data
########################
my $fq1gz = gzopen($fq1,"rb") or die "error open $fq1";
my $fq2gz = gzopen($fq2,"rb") or die "error open $fq2";
open OUT1,'>',"$out\_1.fq" or die "can't open $out\_1.fq\nDied ";
open OUT2,'>',"$out\_2.fq" or die "can't open $out\_2.fq\nDied ";
my $num = 0;
my $clean = 0;
my $outline_1="";
my $outline_2="";
my %Seq=();

while( 1 ){
	my (@q1,@q2);
	for (my $i=0;$i<4 ;$i++) {
		$fq1gz->gzreadline($q1[$i]);
		$fq2gz->gzreadline($q2[$i]);

		if (defined $q1[$i] && defined $q2[$i]){
			chomp ($q1[$i],$q2[$i]);
		}else{
			last;
		}
	}
	unless (defined $q1[3] && defined$q2[3]){
		last;### 当 fq 读完时 退出循环
	}
	$num++;

	if ( exists $Seq{"$q1[1]$q2[1]"} ){
		$Seq{"$q1[1]$q2[1]"} ++;
	#}elsif ( exists $Seq{"$q2[1]$q1[1]"} ){### 考虑到 fq1 fq2 本身没有顺序### 考虑 PCR 后的测序,如果是由 PCR 造成的duplication ,那么测得的 read1 reads2 的顺序应当一致,
		#$Seq{"$q2[1]$q1[1]"} ++;
	}else{
		$outline_1.=join("\n",@q1)."\n";
		$outline_2.=join("\n",@q2)."\n";
		$clean ++;;
		$Seq{"$q1[1]$q2[1]"} =1;
		if ($clean%40000==0) {### 减少输出次数，降低 io
			print OUT1 $outline_1;
			print OUT2 $outline_2;
			$outline_1="";
			$outline_2="";
		}
	}
}
if ($clean%40000!=0) {
	print OUT1 $outline_1;
	print OUT2 $outline_2;
}
$fq1gz->gzclose();
$fq2gz->gzclose();
close OUT1;
close OUT2;

my $Duplicate = 0;
open DIRTY,'>',$out.'.dup.list' or die "can't open the dup list $out.dup.list ";
foreach my $reads ( keys %Seq ){
	print DIRTY $reads."\t".$Seq{$reads}."\n";
	if ($Seq{$reads} > 1) {
		$Duplicate +=$Seq{$reads};
	}
}
close DIRTY;

print OUT "Total_reads:".($num)."\n"."Duplicate_reads:".($Duplicate)."\n";
print OUT "Clean_reads:".($clean)."\n";
my $clean_per=$clean/$num;
#############################added by dushuo#####
if ($clean_per<$Rate){
	open ERR_ADD,">$out.err";
	print ERR_ADD "\nLOW:the dup output_rate is too low--$clean_per lower than $Rate\n";
}
#######################
print OUT "Clean_rate :".$clean_per."\n";
print OUT "END time: ".time."\n";
close OUT;
