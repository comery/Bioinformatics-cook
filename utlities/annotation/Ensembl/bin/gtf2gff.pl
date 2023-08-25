#!/usr/bin/perl -w
use strict;
die "Usage: <gtf file>\n" unless @ARGV == 1;

my %cdsPos;
if ($ARGV[0] =~ /\.gz$/) {
	open IN, "gunzip -c $ARGV[0] | ";
} else {
	open IN, $ARGV[0];
}
while (<IN>) {
	chomp;
	next if /^#/;
	my @info = split /\t/;
	if ($info[1] !~ "pseudogene" && $info[2] eq "CDS") {
		die unless $info[-1] =~ /protein_id\s+"(\S+)\s?";/;
		my $gene_id = $1;
		($info[3], $info[4]) = sort {$a <=> $b}($info[3], $info[4]);
		push @{$cdsPos{$gene_id}}, [$info[3], $info[4], $info[0], $info[6], $info[7]];
	}
}
close IN;

foreach my $gene (keys %cdsPos) {
	@{$cdsPos{$gene}} = sort {$a->[0] <=> $b->[0]} @{$cdsPos{$gene}};
	my ($gene_bg, $gene_ed, $chr, $strand) = ($cdsPos{$gene}->[0]->[0], $cdsPos{$gene}->[-1]->[1], $cdsPos{$gene}->[0]->[2], $cdsPos{$gene}->[0]->[3]);
	my $chr_len = length($chr);
	if ($chr_len <= 2) {
		$chr = "chr$chr";	
	}
	print "$chr\tprotein_coding\tmRNA\t$gene_bg\t$gene_ed\t.\t$strand\t.\tID=$gene;\n";
	foreach my $p (@{$cdsPos{$gene}}) {
		my ($cds_bg, $cds_ed, $chr, $strand, $phase) = @$p;
		my $chr_len = length($chr);
		if ($chr_len <= 2) {
			$chr = "chr$chr";	
		}
		print "$chr\tprotein_coding\tCDS\t$cds_bg\t$cds_ed\t.\t$strand\t$phase\tParent=$gene;\n";
	}
}
