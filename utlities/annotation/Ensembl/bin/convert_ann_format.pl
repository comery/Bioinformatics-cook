#!/usr/bin/perl -w
use strict;

die "Usage: <ann_file> <gff_file>!\n" unless @ARGV == 2;
my $ann_file = shift;
my $gff_file = shift;

my (%hash1, %hash2);
open AN, "gunzip -c $ann_file |";
while (<AN>) {
	chomp;
	next if (/^Ensembl/);
	my @info = (split /\t/);
	$hash1{$info[1]} = $info[0];
	my $go;
	if (@info > 4) {
		$go = $info[4];
		$hash2{$info[0]}{$info[4]} ++;
	}
}
close AN;

print "Protein ID\tGene ID\tGO Term Accession\n";
open GF, $gff_file;
while (<GF>) {
	chomp;
	my $name = (split /\s+/)[-1];
	my $id;
	if ($name =~ /ID=([^;]+);/) {
		$id = $1;
		$id = (split /_/, $id)[-1];
		my $ge = $hash1{$id};
		my @arr;
		die "$ge" unless $hash2{$ge};
		foreach my $g (keys %{$hash2{$ge}}) {
			push @arr, $g;
		}
		my $gos = join(";", @arr);
		print "$id\t$ge\t$gos\n";
	}
}
close GF;


