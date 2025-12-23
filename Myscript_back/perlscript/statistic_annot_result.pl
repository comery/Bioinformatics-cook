#!/usr/bin/perl
##################################################################
#This script is for chloroplast/mitochondria annotation pipoline to 
#statistic the annot result
###################################################################
use strict;
use warnings;
die "Usage:\n\tperl $0 <*.solar.genewise.gff.cds.position.cds> <*.rename.fa> |<STDOUT>" unless (@ARGV == 2 );
open IN,shift;		# *.solar.genewise.gff.cds.position.cds
open LEN,shift;		#*.rename.fa
open OUT,">statistic_report.txt";
open OUT1,">brife_report.txt";
my %hash;
$/=">";<IN>;$/="\n";
while(my $id = <IN>){
	chomp $id;
	my @arr = split /\s+/,$id;
	my $ref = $arr[0];
	my $locus = $1 if ($arr[2] =~ /locus=(\w+):+?\d+/);
   # print "$locus\t$genename\n";
	$/=">";
    my $seq=<IN>;
    chomp $seq;
    $/="\n";
	$hash{$locus}{$ref} = $seq;
}
my %length;
$/ = ">";<LEN>;$/ = "\n";
while (my $id = <LEN>) {
	chomp $id;
	$/ = ">";
	my $seq = <LEN>;
	chomp $seq;
	my $len = length $seq;
	$length{$id} = $len;
	$/ = "\n";
}

foreach my $key1 (keys %hash){
	my $hash2 = $hash{$key1};
	my $number = keys %$hash2;
	print OUT ">$key1\t$number\tlen=$length{$key1}\n";
	print OUT1 ">$key1\t$number\tlen=$length{$key1}\n";
	foreach my $key2 ( keys %$hash2){
		print OUT "\t$key2\n";
		my $gene_name = &gene($key2);
		print OUT1 "$gene_name\t";
	}
	print OUT1 "\n";
}
close IN;
close OUT;
close OUT1;;



sub gene {
	my $str = shift;
	my $genename = (split /\_/,$str)[2];
	return $genename;
}
