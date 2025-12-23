#!/usr/bin/perl
=head1 Description
	This script is for chloroplast/mitochondria annotation pipoline to 
	statistic the annot result. It will output four basic files :
	statistic_report.txt brife_report.txt sorted_annot_result.txt log
	if results can meet the damand, there will be a more file which has
	best matched  scaffold.

=head1 Author :
			YangChentao yangchentao@genomics.cn 2015/04/21

=head1 options
	--cds 	<str>		*.solar.genewise.gff.cds.position.cds
	--fa	<str>		*.fa
	--type <str>		mt(mitochondria) or cp(chloroplast)
	--help 				display this help information

=head1 Usage 
	perl statistic_and_sort_annot_result.pl  -cds <*.solar.genewise.gff.cds.position.cds> -fa <*..fa> -type (mt|cp)

=head1 attention
	So, if you still ask me what is hell of this? I appreciate you take 
	some time to know a little about pipeline for annotating mitochondria
	genome writen by YY(liyiyuan),and then take a look of this script.

=cut

use strict;
use warnings;
use Getopt::Long;

my ($Help,$type,$annot_result,$fasta);
GetOptions(
			"help" => \$Help,
			"cds:s" => \$annot_result,
			"fa:s" => \$fasta,
			"type:s" => \$type
			);
#die `pod2text $0` if (@ARGV==0);
die `pod2text $0` if ($Help || !@ARGV==2 || !$annot_result);
die "Please give a genome type mitochondria or chloroplast. use -type (mt|cp) " unless ($type);
open IN, $annot_result;
open LEN, $fasta;
open OUT,">statistic_report.txt";
open OUT1,">brife_report.txt";
open OUT2,">sorted_annot_result.txt";
open LOG,">log";

###########Read annotated result############
my (%title,%seq,$locus,$start,$end,$middle,%abun);
$/=">";<IN>;$/="\n";
while(my $id = <IN>){
	chomp $id;
	my @arr = split /\s+/,$id;
	my $ref = $arr[0];
	$locus = $1, $start = $2,$end = $3 if ($arr[2] =~ /locus=(\w+):(\d+):(\d+)/);
	$middle = ($start + $end)/2;
	$/=">";
    my $seq=<IN>;
    chomp $seq;
    $/="\n";
	if (defined $title{$locus}{$middle}) {
		$middle = $middle +0.5;
		$title{$locus}{$middle} = $id;
		$seq{$locus}{$middle} = $seq;
	}else {
		$title{$locus}{$middle} = $id;
		$seq{$locus}{$middle} = $seq;
	}
}

#############Read assembly result###########
my (%length,%scaff,@arr,$abundance);
$/ = ">";<LEN>;$/ = "\n";
while (my $id = <LEN>) {
	chomp $id;
	@arr= split /\s+/,$id;
	my $nid = $arr[0];
	if ($nid =~ /scaffold/) {
		$abundance = $arr[2] ;
	}else {
		 $abundance = $arr[1];
	}
#	print "$abundance\n";
	$abun{$nid} = $abundance;
	$/ = ">";
	my $seq = <LEN>;
	chomp $seq;
	my $len = length $seq;
	$length{$nid} = $len;
	$scaff{$nid} = $seq;
	$/ = "\n";
}

#########Set standrad #####################
my (@aa,@bb,$check_aa_num,$check_len_num,$head);
if ($type eq "mt") {
	$check_aa_num = 10;
	$check_len_num = 12000;
	$head = "Mitochondria";
}else {
	$check_aa_num = 50;
	$check_len_num = 100000;
	$head = "Chloroplast";
}
##print head to brife_report.txt#
print OUT1 "Scaff_ID\tAbundance\tAnnoted_Num\tLen\tGenename\n";

my  @sorted_locus = sort(keys %title);
my $annotated_count = @sorted_locus;
my $most_pro = 0;my $most_pro_id; 
my $longest_len = 0;my $longest_len_id;
foreach my $key1 (@sorted_locus){
	my $hash2 = $title{$key1};
	my @sorted_middle = sort{$a <=> $b} keys %$hash2;
	my $number = keys %$hash2;
	print OUT ">$key1\t$abun{$key1}\t$number\tlen=$length{$key1}\n";
	print OUT1 "$key1\t$abun{$key1}\t$number\t$length{$key1}\t";
		push @aa,$key1 	if ($number >= $check_aa_num) ;			## you can modify here to meet your expectation(13 if mitochondria)
		if ($number > $most_pro) {
			$most_pro = $number;
			$most_pro_id = $key1;
		}
		push @bb,$key1	if ($length{$key1} >= $check_len_num);	## maybe 12000 more or less if mitochondria
		if ($length{$key1} >$longest_len) {
			$longest_len = $length{$key1};
			$longest_len_id = $key1;
		}
	foreach my $key2 (@sorted_middle){
		my $genename = &gene($title{$key1}{$key2});
		print OUT "\t$title{$key1}{$key2}\n";
		print OUT1 "$genename\t";
		print OUT2 ">$title{$key1}{$key2}\n$seq{$key1}{$key2}\n";
	}
	print OUT1 "\n";

}
##

####### make a log file#########

my @top;
my $pro = @aa; my $lens = @bb;
print LOG "[$head]\n";
if ($pro == 0 && $lens == 0) {
	print LOG "The number of  scaffold annotated more than $check_aa_num proteins:\t0\nThe number of scaffold's length longer than $check_len_num:\t0";
}elsif ($pro == 0 && $lens >0) {
	@top = @bb;
	print LOG "The number of  scaffold annotated more than $check_aa_num proteins:\t0\nThe number of scaffold's length longer than $check_len_num:\t$lens"
}elsif ($lens == 0 && $pro >0) {
	@top = @aa;
	print LOG "The number of  scaffold annotated more than $check_aa_num proteins:\t$pro\nThe number of scaffold's length longer than $check_len_num:\t0";
}else {
	push @aa,@bb;
	@top = @aa;
	print LOG "The number of  scaffold annotated more than $check_aa_num proteins:\t$pro\nThe number of scaffold's length longer than $check_len_num:\t$lens";
}
##
print LOG "\nThe number of annotated scaffolds:\t$annotated_count\n";
print LOG "The most annotated scaffold's id:\t$most_pro_id\t$most_pro\n";
print LOG "The longest annotated scaffold's id:\t$longest_len_id\t$longest_len\n";

my %count;
my @uniq_top = grep { ++$count{ $_ } < 2; } @top;
if (@uniq_top > 0){
	open MAX, ">max_length_proNum.fa";
	foreach (@uniq_top){print MAX ">$_\n$scaff{$_}";}
	close MAX;
}

close IN;
close OUT;
close OUT1;
close OUT2;
close LOG;

sub gene {
	my $str = shift;
	my $ref =  (split /\s+/,$str)[0];
	my $genename = (split /\_/,$ref)[3];
	return $genename;
}
