#!/usr/bin/perl -w
use strict;
die "Usage : perl $0 <*.mito.fa>" unless (@ARGV == 1);

open IN,shift;

$/=">";<IN>;$/="\n";
while (<IN>){
	chomp;
	my $id = $_;
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/\s//g;
	my $head = substr($seq,0,100);
	my $tail = substr($seq,-100);
	print "head :$head \n";
	print "tail :$tail \n";
	$/="\n";
}

