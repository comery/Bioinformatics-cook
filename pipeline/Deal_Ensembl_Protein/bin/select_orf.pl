#!/usr/bin/perl -w
use strict;
die "Usage: <best.pep> <all.orf>\n" unless @ARGV == 2;

my %list;
open IN, $ARGV[0];
while (<IN>) {
	next unless /^>/;
	my $id = (split /\s+/)[0];
	$id =~ s/^>//;
	$list{$id} ++;
}
close IN;

open IN, $ARGV[1];
while (<IN>) {
	my $id = (split /\s+/)[0];
	print if $list{$id};
}
close IN;
