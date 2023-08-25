#!/usr/bin/perl -w
use strict;
###########
# date 2010-01-03

die "Usage: perl $0 <sort.gff> \n" if @ARGV <1;

my $file = shift;
open IN,"<$file" or die "$!\n";

my ($start,$end);
my $name = "";
while (<IN>) {
	chomp;
	next if (/#/);
	my @t = split /\s+/;
	if ($t[0] eq $name) {
		if ($t[2] <= $end + 1) {
			$end = ($t[3] > $end) ? $t[3] : $end;
		}elsif ($t[2] > $end) {
			print "$name\t$start\t$end\n";
			$name = $t[0];
			$start = $t[2];
			$end = $t[3];
		}	
	}elsif ($t[0] ne $name) {
		if ($name) {
			print "$name\t$start\t$end\n";
		}
		$name = $t[0];
		$start = $t[2];
		$end = $t[3];
	}
}
print "$name\t$start\t$end\n";
close IN;
