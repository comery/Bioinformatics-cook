#!/usr/bin/perl -w
use strict;
die "perl $0 <*.fa>" unless (@ARGV >0) ;
my $genome=shift;
open (IN,$genome) or die $!;
$/=">";<IN>;$/="\n";
while (<IN>) {
	chomp;
	my $id = $_;
	print ">$id\n";
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq =~ s/\n//g;
	my $length=length $seq;
	for (my $i=0; $i<$length; $i+=60) {
		my $part=substr($seq,$i,60);
		print "$part\n";
	}
	$/="\n";
}
close IN;
	
