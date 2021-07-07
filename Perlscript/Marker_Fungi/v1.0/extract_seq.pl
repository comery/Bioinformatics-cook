#!/usr/bin/perl
use strict;
use warnings;
open IN,shift or die "$!";
open IN1,shift or die "$!";
my @tmp;
while (<IN1>) {
chomp;
push @tmp,(split /\s+/,$_)[0];
#print "@tmp";
}

my %hash;
$/=">";<IN>;$/="\n";
while(<IN>){
	chomp;
	my $id=(split,$_)[0];
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\n//g;
	$/="\n";
	$hash{$id}=$seq;
#	print "$id\n$hash{$id}";
}

foreach my $k(@tmp){
	if ($k=~/$id/){
	print ">$k\n$hash{$id}\n";
	}
}
close IN;
close IN1;
