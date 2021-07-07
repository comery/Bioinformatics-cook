#!/usr/bin/perl -w
use strict;
die "Usage: <in file> <sp> <file format(gff or fa)>\n" unless @ARGV == 3;
my $in_file = shift;
my $sp = shift;
my $format = shift;


if ($format eq "gff") {
	&gff_parser($in_file, $sp);
} elsif ($format eq "fa") {
	&fasta_parser($in_file, $sp);
}


sub gff_parser {
my $in_file = shift;
my $sp = shift;
open IN, $in_file;
while (<IN>) {
	my @info = split /\t/;
	my $gene_id;
	$gene_id = $1 if $info[8] =~ /ID=(.+?);/;
	$gene_id = $1 if $info[8] =~ /Parent=(.+?);/;
	my $new_id = "${sp}_$gene_id";
	die "$gene_id\terror!" unless $info[8] =~ s/$gene_id/$new_id/;
	my $out = join "\t", @info;
	print $out;
}
close IN;
}

sub fasta_parser {
my $in_file = shift;
my $sp = shift;
open IN, $in_file;
while (<IN>) {
	if (/^>/) {
		my $id = (split /\s+/)[0];
		$id =~ s/^>//;
		my $new_id = "${sp}_$id";
		die "$id\terror!" unless s/$id/$new_id/;
		print;
	} else {
		print;
	}
}
close IN;
}
