#!/usr/bin/perl
use strict;
use warnings;
open IN,shift or die "$!";
$/="//";<IN>;$/="\n";
print "# STOCKHOLM 1.0\n\n";
while (<IN>){
	chomp;
	if (! $_=~ /^$/) {
		my @aa=split;
		my $id=$aa[0];
		my $seq=join("",@aa[1..$#aa]);
		$seq=~s/\s//g;
		print "$id\t$seq\n";
	}else {
		print "$_\n";
	}
#	print "$_\n";
	}
print "//\n";
close IN;
