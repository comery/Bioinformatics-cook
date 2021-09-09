#!/usr/bin/perl -w
use strict;
die "Usage: <genome fasta>\n" unless @ARGV == 1;

if ($ARGV[0] =~ /\.gz$/) {
	open IN, "gunzip -c $ARGV[0] | ";
} else {
	open IN, $ARGV[0];
}
while (<IN>) {
	if (/^>(\S+)/) {
		my $id = $1;
		my $len = length($id);
		if ($len <= 2) {
			$id = "chr$id";
		}
		print ">$id\n";
	} else {
		print;
	}
}
close IN;
