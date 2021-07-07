#!/usr/bin/perl
use strict;
use warnings;
open IN,shift or die "$!";
$/="//";<IN>;$/="\n";
print "# STOCKHOLM 1.0\n\n";
while (<IN>){
	chomp;
	if (/\S/) {
	#	print "$_\n";
		my @aa=split /\s+/,$_;
		my $id = shift @aa;
		my $seq=join("",@aa[0..$#aa]);
	#	printf "%-50s%-20s\n",$id,$seq;
		print "$id\t$seq\n";
	}else {
		print "$_\n";
	}
}
print "//\n";
close IN;
