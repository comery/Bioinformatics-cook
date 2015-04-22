#!usr/bin/perl -w 
use strict;
open IN,shift;
my ($name1,$name2);
while (<IN>) {
	chomp;
	my $len = $_;
	#print $len;
	#$name1 = $1 if ($len =~ m/\_(\w+?)\_/);
	$name2 = $1 if ($len =~ m/(links\S+\.txt$)/);
	print "$name2";
	}
close IN;
