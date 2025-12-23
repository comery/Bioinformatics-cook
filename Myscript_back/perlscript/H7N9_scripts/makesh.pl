#!/usr/bin/perl -w 
use strict;
open IN,shift;
while (<IN>) {
	chomp;
	my @a = split;
	my $inputfile = $a[0];
	my $start = $a[2];
	my $end = $a[3];
	my $gene = $a[1];
	print "perl /ifs4/NGB_ENV/USER/yangchentao/project/new_five_sample_20150508/scripts/syn_and_nsyn.pl $inputfile -start $start -end $end   -pro $gene\n";
}
close IN;
