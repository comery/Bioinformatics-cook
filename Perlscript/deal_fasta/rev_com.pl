#!/usr/bin/perl -w 
use strict;
die "Usage:perl $0 <*.fa>" unless (@ARGV > 0);
open IN, shift;
$/=">";<IN>;$/="\n";
while (my $id = <IN>) {
	chomp $id;
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/\n//g;
	$seq =~ tr/TAGCatcg/ATCGtagc/;
	$seq= reverse($seq);
	print ">$id\n";
	&wrapseq($seq);
	$/="\n";
}

sub wrapseq{
	my $seq = shift;
	my $number = 80;
	my $length = length($seq);
	for (my $i=0; $i<$length; $i+=$number) {
		my $part=substr($seq,$i,$number);
		print "$part\n";
	}
}
