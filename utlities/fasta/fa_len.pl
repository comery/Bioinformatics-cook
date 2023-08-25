#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

die "usage: perl $0 contig.fasta [contig1.fasta] [contig2.fa]\n" unless (@ARGV>0);


foreach my $f(@ARGV) {
	if ($f =~ /gz$/){
		open IN,"gzip -dc $f |" || die "can't open the input file,$!";
	}else{
		open IN,$f or die "$!\n";
	}

	open OUT,">$f.length" or die "$!\n";
	my($head,$seq);
	my @line;
	my $line;
	$/=">";<IN>;$/="\n";
	my (%length, %seq_length);

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
	
		$length{$head} = $len;
		$seq_length{$head} = $seq_len;
		$/="\n";
	}

	my @sorted = sort{$length{$b} <=> $length{$a}} keys %length;
	foreach (@sorted){
		print OUT "$_\t$length{$_}\t$seq_length{$_}\n";
	}
	my $count = keys %length;

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(values %length);
	print "-------".$f."-------\n";
	print "total scaffold number:\t";
	print "$count\n";
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
	close IN;
	close OUT;
}
