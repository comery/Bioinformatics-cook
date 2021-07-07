#!/usr/bin/perl -w
use strict;
use Getopt::Long;
die "Usage: perl $0 <fasta> 150 " unless (@ARGV == 2) ;

open IN,shift;
my $l = shift;

$/=">";<IN>;$/="\n";
while (my $id = <IN>) {
	chomp $id;
	$/= ">";
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/\n//g;
	my $len = length $seq;
	if ($len >= $l) {
		1 while $seq =~ s/(\w{60})(\w+)$/$1\n$2/;
		print ">$id\n$seq\n";
	}
	$/="\n";
	
}

close IN;
