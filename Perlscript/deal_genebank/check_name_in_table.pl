#!/usr/bin/perl -w 
use strict;
my @names = @ARGV;
die "Give me a mitochondrial gene name to check!" unless (@names > 0);
my %hash;
open IN ,"/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/perlscript/deal_genebank/MT_Pname.hash.txt" or die "$!";
while (<IN>) {
	chomp;
	next if (/^#/);
	my @a = split /\|/,$_;
	my $gene = $a[0];
	my @iso = split /,/,$a[1];
#	print "@iso\n";
	foreach my $n(@iso){
		$hash{$n} = $gene if ($n);
	}
}

foreach my $key(@names) {
	chomp $key;
	if (exists $hash{$key}) {
		print "$key in table!\t$key => $hash{$key}\n";
	}else {
		print "$key not in table!\n";
	}
}

close IN;
