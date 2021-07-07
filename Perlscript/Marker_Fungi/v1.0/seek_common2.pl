#!usr/bin/perl -w 
use strict;
my (@a,$id,$pro,%hash,%count);
die "Usage:\t\nperl $0 <all.id> <all.fa> <otudir>|<STDOUT>" unless (@ARGV == 3);
open IN,$ARGV[0];#all_best_to_best.id
open FA,$ARGV[1];
my $outdir = $ARGV[2];
`mkdir $outdir` unless (-d "./$outdir");
open OUT, ">orthology_distribution1.xls";
open OUT1, ">common.id1.list";
while (my $str=<IN>) {
	@a = split /\s+/,$str;
	$id = $a[0];
	$pro = $a[1];
	if (defined $hash{$pro}) {
		$hash{$pro} .= "\t$id";
		$count{$pro} ++;
	}else{
		$hash{$pro} = $id;
		$count{$pro} =1;
	}
}
my %cds;
$/ = ">";<FA>;$/ = "\n";
while (<FA>) {
	chomp ;
	my $title = $_;
	$/ = ">";
	my $seq = <FA>;
	chomp $seq;
	$cds{$title} = $seq;
	$/ = "\n";
}

foreach my $key (keys %count) {
	print OUT "$key\t$hash{$key}\n";
	if ($count{$key} == 53) {
		print OUT1 "$key\n" ;
		my @b = split /\t/,$hash{$key};
		`mkdir $outdir/$key` unless (-d "$outdir/$key");
		open CD ,">$outdir/$key/$key.fa";
		foreach (@b) {print CD ">$_\n$cds{$_}";}
		close CD;
	}

}
close IN;
close OUT;
close OUT1;
