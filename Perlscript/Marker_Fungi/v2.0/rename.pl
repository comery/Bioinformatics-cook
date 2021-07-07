#!/usr/bin/perl
use strict;
use warnings;
open IN,shift;
open OUT,">id.list";
my %hash;
$/=">";<IN>;$/="\n";
while(<IN>){
	my $id=$_;
#	chomp $id;
	my $nid=(split /\s+/,$id)[0];
		chomp $nid;
        $/=">";
        my $seq=<IN>;
        chomp $seq;
        $/="\n";
	$hash{$nid}=$seq;
#	print "$seq\n";
	print ">$nid\n$hash{$nid}";
	print OUT "$nid\n";
}
close IN;
close OUT;
