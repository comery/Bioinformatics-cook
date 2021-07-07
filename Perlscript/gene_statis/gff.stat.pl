#!/usr/bin/perl -w
use strict;
die "Usage: perl <*.gff> <type: e.g. Denovo>\n" unless @ARGV == 2;
my $type = $ARGV[1];

my ($total,$Single_exon_gene,$gene_len,$mRNA_len,$exons_num,$exons_len,$intron_len,$intron_num);
my %hash;
open (IN,$ARGV[0]) or die $!;
while (<IN>) {
	chomp;
	my @info = split /\s+/;
	next if ($info[2] =~ /UTR/i || $info[2] =~ /mRNA/i);
	die "Wrong!!!\t$_\n" unless ($info[2] =~ /cds/i);
	my $id = (split /=/,(split /;/,$info[8])[0])[1];
	push @{$hash{$id}},[$info[3],$info[4]];
}
close IN;

foreach my $geneid (sort keys %hash) {
	$total ++;

	if (@{$hash{$geneid}} == 1) {
		$Single_exon_gene ++;
	}

	@{$hash{$geneid}} = sort {$a->[0] <=> $b->[0]} @{$hash{$geneid}};

	foreach my $ref (@{$hash{$geneid}}) {
		my ($st,$ed) = @{$ref};
		my $cds = $ed-$st+1;
		$mRNA_len += $cds;
		$gene_len += $cds;
		$exons_num ++;
	}

	for (my $i = 0; $i < @{$hash{$geneid}}-1; $i ++) {
		my ($bg1, $ed1) = @{$hash{$geneid}->[$i]};
		my ($bg2, $ed2) = @{$hash{$geneid}->[$i+1]};
		my ($intron_bg, $intron_ed) = (sort {$a <=> $b}($bg1, $ed1, $bg2, $ed2))[1,2];
		$intron_bg ++;
		$intron_ed --;
		$intron_len += $intron_ed-$intron_bg+1;
		$gene_len += $intron_ed-$intron_bg+1;
		$intron_num ++;
	}
}

my $ave_gene = sprintf "%.0f",$gene_len/$total;
my $ave_mRNA = sprintf "%.0f",$mRNA_len/$total;
my $ave_exon =  sprintf "%.2f",$exons_num/$total;
my $ave_exon_len = sprintf "%.0f",$mRNA_len/$exons_num;
my $ave_intron = sprintf "%.0f",$intron_len/$intron_num;
my $line = join "\t",$type,$total,"unknown",$Single_exon_gene,$ave_gene,$ave_mRNA,$ave_exon,$ave_exon_len,$ave_intron;
print "$line\n";
