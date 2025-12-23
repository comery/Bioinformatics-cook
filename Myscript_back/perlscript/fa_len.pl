#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

die "usage: perl $0 [ contig file ]\n" unless (@ARGV==1);
open IN,$ARGV[0] or die "$!\n";
open OUT,">$ARGV[0].length" or die "$!\n";
my($head,$seq);
my @line;
my $line;
$/=">";<IN>;$/="\n";
my (@length, @seq_length);
#my ($max,$min);
#$max = 0;
#$min = 100000;

while ($line=<IN>) {
	@line=split(/\s+/,$line);
        my $head=$line[0];
        chomp $head;        
        $/=">";
        my $seq=<IN>;
        chomp $seq;
	$seq=~s/\s+//g;
	$seq=~s/>//g;  
    	my $len=length($seq);
	$seq =~ s/N//g;#please mofity this step if you want check the protein seq length!
	my $seq_len = length($seq);
#
	
	push @length,$len;
	push @seq_length,$seq_len;
     	print OUT "$head\t$len\t$seq_len\n";
        $/="\n";

#	if($len<$min){
#		$min = $len;
#	}elsif($len>$max){
#		$max = $len;
#	}
}
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@length);
print "minimal value\t";
print $stat->quantile(0); # => 1
print "\nlower quartile\t";
print $stat->quantile(1); # => 3.25
print "\nmedian\t";
print $stat->quantile(2); # => 5.5
print "\nupper quartile\t";
print $stat->quantile(3); # => 7.75
print "\nmaximal value\t";
print $stat->quantile(4); # => 10
print "\n";
