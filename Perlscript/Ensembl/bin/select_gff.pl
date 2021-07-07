#!/usr/bin/perl -w
use strict;
die "Usage: <gene list> <gff files>\n" unless @ARGV >= 2;
my $list = shift;
my %gff;
foreach my $in_file (@ARGV) {
	&getGff($in_file, \%gff);
}

open IN, $list;
while (<IN>) {
	next if /^#/;
	my @info = split /\s+/;
	my $id = $info[0];
	if ($gff{$id}) {
		print $gff{$id};
	} else {
		print STDERR "$id is miss in gff file\n";
	}
}
close IN;

sub getGff {
	my ($in_file, $ref) = @_;
	open IN, $in_file;
	while (<IN>) {
		my @info = split /\s+/;
		my $id;
		if ($info[-1] =~ /^ID=(\S+?);/) {
			$id = $1;
		} elsif ($info[-1] =~ /^Parent=(\S+?);/) {
			$id = $1;
		} else {
			die;
		}
		$ref->{$id} .= $_;
	}
	close IN;
}
