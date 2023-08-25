#!/usr/bin/perl -w
use strict;

die "perl $0 <in.fa> <out.fa> " if @ARGV!=2;

open(IN, "<", $ARGV[0]) or die "$!";
open(OUT, ">", $ARGV[1]) or die "$!";
my ($ti, $seq);
$/ = "\>";
<IN>;
$/ = "\n";
while($ti=<IN>){
	print OUT ">$ti";
	
	$/ = "\>";
	chomp($seq=<IN>);
	$/ = "\n";

	$seq =~ s/\s+//g;
	$seq =~ tr/atgcn/NNNNN/;
	1 while $seq =~ s/(\w{60})(\w+)$/$1\n$2/;
	print OUT "$seq\n";
}
close IN;
close OUT;
