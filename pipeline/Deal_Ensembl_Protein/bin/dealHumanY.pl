#!/usr/bin/perl -w
use strict;
use lib '/share/project002/liqiye/bin/module/personal';
use Fasta;
die "Usage: <human.fa>\n" unless @ARGV == 1;

if ($ARGV[0] =~ /\.gz$/) {
	open IN, "gunzip -c $ARGV[0] | ";
} else {
	open IN, $ARGV[0];
}
$/ = ">";
<IN>;
while (<IN>) {
	chomp;
	/(.+)\n/;
	my $id = (split /\s+/, $1)[0];
	if ($id eq "chrY") {
		s/.+\n//;
		s/\s+|>//g;
		my $len = length ($_);
		my $n = ($_ =~ s/N/N/ig);
		my $n_ratio = $n/$len;
		next if $n_ratio > 0.8;
		print ">$id\n";
		$_ = sequence_out($_, 60);
		s/\s+$//;
		print "$_\n";
	} else {
		print ">$_";
	}
}
$/ = "\n";
close IN;
