#!/usr/bin/perl
use strict;
use warnings;
#open IN1,shift or die "$!";#a fam_id of com_all_best_to_best.id
open IN,shift or die "$!";#all.fas
my @tmp;
while (<>) {
chomp;
push @tmp,$_;
}

my %hash;
$/=">";<IN>;$/="\n";
while(<IN>){
	chomp;
	my $id=$_;
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$/="\n";
	$hash{$id}=$seq;

foreach my $k(@tmp){
	if ($id=~/$k/){
	print ">$id\n$hash{$id}";
	}
}
}
close IN;
#close IN1;
