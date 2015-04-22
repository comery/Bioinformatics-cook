#!usr/bin/perl -w 
use strict;
open IN,shift;
my $line = <IN>;
$line =~ s/\s+//g;
my $replaced = substr ($line,22);
#print "$replaced\n";
my $len = length($replaced);
my $i;
for ($i=0;$i <= $len-3;$i+=2) {
	my $tmp=substr ($replaced,$i,2);
	my @singal = split //,$tmp;
	print "@singal\n";
	}
close IN;
