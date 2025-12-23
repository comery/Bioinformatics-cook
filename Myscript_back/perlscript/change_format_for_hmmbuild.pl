#!/usr/bin/perl
use strict;
use warnings;
open IN,shift or die "$!";
$/="//";<IN>;$/="\n";
print "# STOCKHOLM 1.0\n\n";
while (<IN>){
	chomp;
	my @aa=split;
	my $id = $aa[0];
	#print "$id\n";
	my $seq=join("",@aa[1..$#aa]);
	
#	$seq=~s/\s//g;
#	print "$id\t$seq\n";
	printf "%-50s%-20s\n",$id,$seq;
	}
print "//\n";
close IN;
