#!/usr/bin/perl 
=head1	description

this program is to praise the rRNA seq from transcription data accroding the information
from annotation and fpkm.

=head1	Author & Vesion
	Author: Chentao YANG, yangchentao@genomics.org.cn
	Version: 1.0,  Date: 2014-12-20
	
=head1	options
	
	--fa	<string> 	fasta file
	--b	<string>	blast.out
	--fp	<string>	fpkm.xls
	--o 	<string>	output dir
	--top	<number>	the number of best hits
	--help	print this help

=head1 Usage
	
	perl praise_and_output.pl -fa Unigene.fa -b blast.out -fp fpkm.xls -o out -top 5

=cut

#########################################################################################

use strict;
use Getopt::Long;
my ($Blast,$Outdir,$Scaff,$fpkm,$Help,$top);
GetOptions(
	"fa:s"=> \$Scaff,
	"b:s"=> \$Blast,
	"o:s"=> \$Outdir,
	"fp:s"=> \$fpkm,
	"top:n"=> \$top,
	"help"=> \$Help
	);
$Outdir ||= "./";
$top ||= 5;
#`mkdir $Outdir` if ($Outdir ne "./");
die `pod2text $0` if (@ARGV == 0||$Help);
open FA,"$Scaff";#scaffold assembly
open BLA,"$Blast";#blast.out 
open FPKM,"$fpkm";#fpkm.xls
open OUT1,">$Outdir/5S.fa";
open OUT2,">$Outdir/5.8S.fa";
open OUT3,">$Outdir/18S.fa";
open OUT4,">$Outdir/28S.fa";
open OUT5,">$Outdir/ITS.fa";
open OUT6,">$Outdir/MIX.fa";
#open OUT7,">$Outdir/rRNA.fa";
$/=">";<FA>;$/="\n";
###############################
my %scaff;
while (my $id=<FA>) {
	chomp $id;
	my $nid = (split /\s+/,$id)[0];
	$/=">";
	my $seq=<FA>;
	chomp $seq;
	$seq=~s/\s//g;
#	print "$nid\n$seq\n";
	$scaff{$nid} = $seq;
	$/="\n";
}
close FA;
###############################
my %info;
while (<FPKM>) {
	chomp;
	next if ($_ =~ /geneID/);
	my @array = split /\s+/,$_;
	my $gid = $array[0];###Be careful here,motify the number for real condition
	my $val = $array[3];
#	print "$gid\t$val\n";
	$info{$gid} = $val;
	}
close FPKM;
##############################
my (%hash5S,%hash5_8S,%hash18S,%hash28S,%hashITS,%complex,%len,$check,$id);
while (<BLA>) {
	next if ($_ =~ "Query_id");
	chomp;
	my @array1 = split /\t/,$_;
	$id = $array1[0];
	$check = $array1[15];
	$len{$id} = $array1[1];
	if ($check =~ /\b5S\b/ && !($check =~ /\b5.8S\b/) && !($check =~ /\b18S\b/) && !($check =~ /\b28S\b/) && !($check =~ /\bITS\d?\b/)) {
		if (!exists $hash5S{$id}) {
			$hash5S{$id} = $check;
		}else{
			$hash5S{$id} = $hash5S{$id}."|".$check;
		}
	}elsif ($check =~ /\b5.8S\b/ && !($check =~ /\b5S\b/) && !($check =~ /\b18S\b/) && !($check =~ /\b28S\b/) && !($check =~ /\bITS\d?\b/)) {
		if (!exists $hash5_8S{$id}) {
			$hash5_8S{$id} = $check;
		}else{
			$hash5_8S{$id} = $hash5_8S{$id}."|".$check;
		}
	}elsif ($check =~ /\b18S\b/ && !($check =~ /\b5S\b/) && !($check =~ /\b5.8S\b/) && !($check =~ /\b28S\b/) && !($check =~ /\bITS\d?\b/)) {
		if (!exists $hash18S{$id}) {
			$hash18S{$id} = $check;
		}else{
			$hash18S{$id} = $hash18S{$id}."|".$check;
		}
	}elsif ($check =~ /\b28S\b/ && !($check =~ /\b5S\b/) && !($check =~ /\b18S\b/) && !($check =~ /\b5.8S\b/) && !($check =~ /\bITS\d?\b/))	{
		if (!exists $hash28S{$id}) {
			$hash28S{$id} = $check;
		}else{
			$hash28S{$id} = $hash28S{$id}."|".$check;
		}
	}elsif ($check =~ /\bITS\d?\b/ && !($check =~ /\b5S\b/) && !($check =~ /\b18S\b/) && !($check =~ /\b28S\b/) && !($check =~ /\b5.8S\b/)) {
		if (!exists $hashITS{$id}) {
			$hashITS{$id} = $check;
		}else{
			$hashITS{$id} = $hashITS{$id}."|".$check;
		}
	}elsif (!($check =~ /\bITS\d?\b/) && !($check =~ /\b5S\b/) && !($check =~ /\b18S\b/) && !($check =~ /\b28S\b/) && !($check =~ /\b5.8S\b/)) {
	next;
	}else{
		if (!exists $complex{$id}) {
			$complex{$id} = $check;
		}else{
			$complex{$id} = $complex{$id}."|".$check;
		
		}
	}
}
my ($key1,$key2,$key3,$key4,$key5,$key6);
my (@array5S,@array5_8S,@array18S,@array28S,@arrayITS);
my ($len1,$len2,$len3,$len4,$len5);
foreach $key6(keys %complex){
	if ($len{$key6}>3000 && $info{$key6}>1000){
		print OUT6 ">$key6\t$info{$key6}\t$complex{$key6}\n$scaff{$key6}\n";
	}
}
$len1 = keys %hash5S ;
$len1 >$top ? (@array5S = &praise (keys %hash5S)) : (@array5S = keys %hash5S);
if (!@array5S == 0) {
	foreach  $key1(@array5S){
		print OUT1 ">$key1\t$info{$key1}\t$hash5S{$key1}\n$scaff{$key1}\n";
	}
}
$len2 = keys %hash5_8S ;
$len2 > $top ? (@array5_8S = &praise (keys %hash5_8S)) : (@array5_8S = keys %hash5_8S);
if (!@array5_8S == 0){
	foreach  $key2(@array5_8S){
		print OUT2 ">$key2\t$info{$key2}\t$hash5_8S{$key2}\n$scaff{$key2}\n";
	}
}
$len3 = keys %hash18S ;
$len3 > $top ? (@array18S = &praise (keys %hash18S)) : (@array18S = keys %hash18S);
if (!@array18S == 0) {
	foreach  $key3(@array18S){
		print OUT3 ">$key3\t$info{$key3}\t$hash18S{$key3}\n$scaff{$key3}\n";
	}
}
$len4 = keys %hash28S ;
$len4 > $top ? (@array28S = &praise (keys %hash28S)) : (@array28S = keys %hash28S);
if (!@array28S == 0) {
	foreach  $key4(@array28S){
		print OUT4 ">$key4\t$info{$key4}\t$hash28S{$key4}\n$scaff{$key4}\n";
	}
}
$len5 = keys %hashITS;
$len5 > $top ? (@arrayITS = &praise (keys %hashITS)) : (@arrayITS = keys %hashITS);
if (!@arrayITS == 0) {
	foreach  $key5(@arrayITS){
		print OUT5 ">$key5\t$info{$key5}\t$hashITS{$key5}\n$scaff{$key5}\n";
	}
}

close BLA;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
#close OUT7;
##################################################
sub praise {
	my ($k,@aa,@bb,@cc,%sub_fpkm,%re);
	foreach $k(@_){
		$sub_fpkm{$k} = $info{$k};
	}
	%re = reverse %sub_fpkm;
	@aa = sort{$b<=>$a} keys %re;
	#print "@aa\n";
	@bb = splice (@aa,0,$top);
	my $i;
	foreach  $i(@bb) {
		push @cc,$re{$i};
	}
	return @cc;
}
