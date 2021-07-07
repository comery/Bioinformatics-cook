#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

die "usage: perl $0 [ contig file ]\n" unless (@ARGV==1);
open IN,$ARGV[0] or die "$!\n";
my($head,$seq);
$/=">";<IN>;$/="\n";
my (%length,@len,%hash);
#my ($max,$min);
#$max = 0;
#$min = 100000;

while ($head=<IN>) {
    chomp $head;        
    $/=">";
    my $seq=<IN>;
    chomp $seq;
	$hash{$head} = $seq;
	$seq=~s/\s+//g;
	$seq=~s/>//g;  
    my $len=length($seq);
#	$seq =~ s/N//g;#please mofity this step if you want check the protein seq length!
	$length{$head} = $len;
#
	push @len,$len;
#	push @seq_length,$seq_len;
 #    	print OUT "$head\t$len\t$seq_len\n";
        $/="\n";

#	if($len<$min){
#		$min = $len;
#	}elsif($len>$max){
#		$max = $len;
#	}
}
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@len);
#print "minimal value\t";
#print $stat->quantile(0); # => 1
#print "\nlower quartile\t";
#print $stat->quantile(1); # => 3.25
#print "\nmedian\t";
my $mean = $stat->quantile(2); # => 5.5
#print "\nupper quartile\t";
#print $stat->quantile(3); # => 7.75
#print "\nmaximal value\t";
#print $stat->quantile(4); # => 10
#print "\n";
#print "$mean\n";
while (my ($key,$val) = each %hash) {
	if (abs($length{$key}-$mean) <= 0.2*$mean ) {
		print ">$key\n$val";
	}
}

close IN;
