#!usr/bin/perl -w
use strict;
die "perl $0 <*.fa> <number of one line>" unless (@ARGV == 2) ;
my $genome=shift;
my $num = shift;
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
	1 while $seq =~ s/(\w{$num})(\w+)$/$1\n$2/;
	print "$seq\n";
	$/="\n";
}
close IN;
	
