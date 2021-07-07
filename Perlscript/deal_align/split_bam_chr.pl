#!/usr/bin/perl -w
use strict;
if (@ARGV < 2){
	print "Usage: perl $0 <bam|whole genome mapping> <outdir>";
	exit()
}
my $file = shift;
my $outdir = shift;
open (IN,"samtools view -h $file |") or die $!;


my (@heads,@group);
my $rg = "\@RG\tID:foo\tPL:illumina\tSM:OFFSP";
my $chr = "C";

while (<IN>){
	chomp;
	my $all = $_;
	if ($all =~ /^@/){
		push @heads,$all;
	}else{
		my @t = split;
		$chr = $t[2] if ($chr eq "C");
		next if ($t[6] ne "=");

		if ($t[2] eq $chr){
			push @group, $all;
		}else{
			my $head = join"\n",@heads;
			my $out = join"\n",@group;
			open (OUT,">$outdir/$chr.sam") or die $!;
			print OUT "$head\n$rg\n$out\n";
			close OUT;
			$chr = $t[2];
			@group = ();
			push @group, $all;
		}
	}
}
close IN;

my $head = join"\n",@heads;
my $out = join"\n",@group;
open (OUT,">$outdir/$chr") or die $!;
print OUT "$head\n$rg\n$out\n";
close OUT;
