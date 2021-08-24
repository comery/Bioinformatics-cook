#!/usr/bin/perl -w
use strict;

use lib '/hwfssz1-tmp/ST_DIVERSITY/PUB/USER/zhouyang/bin/convert/personal/';
use GFF;

die "perl $0 <ref.codon.bed> <snp.txt>" unless @ARGV == 2;
my %snps;
my %var;
open(IN, $ARGV[1]) or die $!;
while (<IN>) {
	my @info = split /\t/;
	my $ref = $info[0];
	my $pos = $info[1];
	my $v = $info[3];
	my $mut = $info[4];
	$snps{$ref}{$pos} = $mut;
	$var{$ref}{$pos} = $v;
}

close IN;

open(IN1, $ARGV[0]) or die $!;
while(<IN1>){
	chomp;
	my $type;
	my ($gene, $ref,$strand, $p1, $p2, $p3, $codon, $prot, $index) = split /\t/;
	my ($pos1, $b1) = split /:/, $p1;
	my ($pos2, $b2) = split /:/, $p2;
	my ($pos3, $b3) = split /:/, $p3;
	if (exists $snps{$ref}{$pos1}){
		my $v = $var{$ref}{$pos1};
		my $mut_codon = $snps{$ref}{$pos1} . $b2 . $b3;
		my $mut_prot = &translate($mut_codon);
		if ($mut_prot eq $prot){
			$type = 'syn';
		}elsif (($prot ne 'U') && ($mut_prot eq 'U')){
			$type = 'stop';
		}else{
			$type = 'non-syn';
		}
		print("$_\t$mut_codon\t$mut_prot\t$type\t$v\n");
	}
	if (exists $snps{$ref}{$pos2}){
		my $v = $var{$ref}{$pos2};
		my $mut_codon = $b1 . $snps{$ref}{$pos2} . $b3;
		my $mut_prot = &translate($mut_codon);
		if ($mut_prot eq $prot){
			$type = 'syn';
		}elsif (($prot ne 'U') && ($mut_prot eq 'U')){
			$type = 'stop';
		}else{
			$type = 'non-syn';
		}
		print("$_\t$mut_codon\t$mut_prot\t$type\t$v\n");
	}
	if (exists $snps{$ref}{$pos3}){
		my $v = $var{$ref}{$pos3};
		my $mut_codon = $b1 . $b2 . $snps{$ref}{$pos3};
		my $mut_prot = &translate($mut_codon);
		if ($mut_prot eq $prot){
			$type = 'syn';
		}elsif (($prot ne 'U') && ($mut_prot eq 'U')){
			$type = 'stop';
		}else{
			$type = 'non-syn';
		}
		print("$_\t$mut_codon\t$mut_prot\t$type\t$v\n");
	}

}

sub translate{

	my %CODE = (
			"standard" =>
				{	
				'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
				'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
				'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
				'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
				'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
				'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
				'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
				'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
				'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
				'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
				'ATG' => 'M',                                                                         # Methionine
				'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
				'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
				'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
				'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
				'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
				'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
				'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
				'TGG' => 'W',                                                                         # Tryptophan
				'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
				'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'                                              # Stop
				}
			## more translate table could be added here in future
			## more translate table could be added here in future
			## more translate table could be added here in future
	);
	my $codon = shift;
	my $table = 'standard';
	my $prot = (exists $CODE{$table}->{$codon}) ? $CODE{$table}->{$codon} : 'X';
	return $prot;
}



