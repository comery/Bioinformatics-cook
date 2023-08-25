#!/usr/bin/perl -w
use strict;
die "Usage:perl $0 <all_cds.posi.list> <All-Unigene.ssr.out.xls>" unless (@ARGV == 2);
open IN , shift;
open SR, shift;
open LOG, ">ssr_stat.log";
my %hash;
while (<IN>) {
	chomp;
	my @a = split /\s+/,$_;
	push @{$hash{$a[0]}},$a[1];
	push @{$hash{$a[0]}},$a[2];
}
my $in = 0;
my $out = 0;
my $on = 0;
my $not = 0;
my $total =0;
while (<SR>) {
	$total++;
	chomp;
	my $str = $_;
	my @b = split /\s+/,$str;
	my $id = $b[0];
	my $old = (split /_/,$id)[0]."_All";
#	print "$old\n";
	if (exists $hash{$old}) {
		my $start = $b[3];
		my $len = $b[4];
		my $end = $start + $len -1;
#		print "$start\t$end\t$len\n";
		if ($start >= $hash{$old}[0] && $end <= $hash{$old}[1]) {
			print "IN\t$str\n";
			$in++;
		}elsif ($end < $hash{$old}[0] || $start > $hash{$old}[1]) {
			print "OUT\t$str\n";
			$out++;
		}else {
			print "repeat on start or end of cds!\t$str\n";
			$on++;
		}
	}else {
		print "NOT_CDS\t$str\n";
		$not++;
	}
}
print LOG "SSR in cds regin:\t$in\n";
print LOG "SSR out of cds regin:\t$out\n";
print LOG "SSR  over cds regin:\t$on\n";
print LOG "SSR not in cds regin:\t$not\n";
print LOG "SSR total number:\t$total";
close IN;
close SR;
		

