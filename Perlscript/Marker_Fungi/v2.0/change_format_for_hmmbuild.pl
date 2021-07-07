#!/usr/bin/perl
use strict;
use warnings;
die "Usage: perl $0 <*.out|msf format>" unless (@ARGV==1);

open IN,shift or die "$!";
$/="//";<IN>;$/="\n";
print "# STOCKHOLM 1.0\n\n";
while (<IN>){
	chomp;
	if (! $_){
		print "\n";
	}else {
		my @aa=split;
		my $id = $aa[0];
		die "The id \"$id\" is too long!" if (length $id >= 50);
		#print "$id\n";
		my $seq=join("",@aa[1..$#aa]);
		printf "%-50s%-20s\n",$id,$seq;
	}
	
}
print "//\n";
close IN;
