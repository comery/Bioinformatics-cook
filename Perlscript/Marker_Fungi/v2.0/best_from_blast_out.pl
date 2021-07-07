#!usr/bin/perl -w 
use  strict;
die "Usage: perl $0 <blastp.out> " unless (@ARGV==1);
open IN,shift;
open OUT,">blastp.best.out";
open OUT1,">blastp.best.hit";
my (%hash,@a,%hit);
while (my $line=<IN>) {
	chomp $line;
	@a=split /\s+/,$line;
	my $gene = $a[0];
	my $fam_id=(split /_/,$a[1])[-1];
#	print "$gene\t$fam_id\n";
	if ( exists $hash{$gene}) {
		next;
	}else {
		$hash{$gene} = $fam_id;
		$hit{$gene} = $line;
	}

}

foreach my $key (keys %hash) {
	print OUT "$key\t$hash{$key}\n";
	print OUT1 "$hit{$key}\n";
}
	
