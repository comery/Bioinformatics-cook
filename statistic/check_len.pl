#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

die "usage: perl $0 [ contig file ]\n" unless (@ARGV==1);
open IN,$ARGV[0] or die "$!\n";
#open OUT,">$ARGV[0].length" or die "$!\n";
my($head,$seq);
my ($line,%hash);
$/=">";<IN>;$/="\n";
my @length;
while ($line=<IN>) {
	$head=(split /\s+/,$line)[0];
    chomp $head;        
    $/=">";
    my $seq=<IN>;
    chomp $seq;
	$seq=~s/\s+//g;
	$seq =~ s/N//g;
	$hash{$head}=$seq;
    my $len=length($seq);
	push @length,$len;
#    print OUT "$head\t$len\t$seq_len\n";
    $/="\n";
}
my	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@length);
my $min=$stat->quantile(0); # => 1
my $max=$stat->quantile(4); # => 10
print "$max\t$min\n";
if (($max-$min)/$max<0.2){
	foreach my $key(keys %hash){
		print ">$key\n$hash{$key}\n";
	}
}
