#!/usr/bin/perl
##################################################################
#This script is for chloroplast/mitochondria annotation pipoline to 
#statistic the annot result
###################################################################
use strict;
use warnings;
open IN,shift;# *.solar.genewise.gff.cds.position.cds
open OUT,">sorted_annot_result.txt";
my (%title,%seq,$locus,$start,$end,$middle);
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
my  @sorted_locus = sort(keys %title);
foreach my $key1 (@sorted_locus){
	my $hash2 = $title{$key1};
	my @sorted_middle = sort{$a <=> $b} keys %$hash2;
	foreach my $key2 ( @sorted_middle){
		print OUT "$title{$key1}{$key2}\n$seq{$key1}{$key2}\n";
	}
}
close IN;
close OUT;
