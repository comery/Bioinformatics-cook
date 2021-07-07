#!/usr/bin/perl -w
use strict;

die "Usage: <ann_file> <gff_file>!\n" unless @ARGV == 2;
my $ann_file = shift;
my $gff_file = shift;

my (%hash1, %hash2, %name);
open AN, "gunzip -c $ann_file |";
while (<AN>) {
	chomp;
	next if (/^Ensembl/);
	my @info = (split /\t/);
	$hash1{$info[1]} = $info[0];
	$info[2] = "-" unless $info[2];
	$info[3] = "-" unless $info[3];
	$name{$info[0]} = [$info[2], $info[3]];
	my $go;
	if (@info > 4) {
		$go = $info[4];
		$hash2{$info[0]}{$info[4]} ++;
	}
}
close AN;

print "1.Ensembl_Protein_ID\t2.Ensembl_Gene_ID\t3.GO_Term_Accession\t4.Associated_Gene_Name\t5.Description\n";
open GF, $gff_file;
while (<GF>) {
	chomp;
	my $name = (split /\t+/)[-1];
	my $id;
	if ($name =~ /ID=([^;]+);/) {
		$id = $1;
		next unless $hash1{$id};
		my $ge = $hash1{$id};
		my @arr;
		if ($hash2{$ge}) {
			foreach my $g (keys %{$hash2{$ge}}) {
				push @arr, $g;
			}
		}
		my $gos = join(";", @arr);
		$gos = "-" unless $gos;
		my $des = join("\t", @{$name{$ge}});
		print "$id\t$ge\t$gos\t$des\n";
	}
}
close GF;


