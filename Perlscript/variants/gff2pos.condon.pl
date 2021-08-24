#!/usr/bin/perl -w
use strict;

use lib '/hwfssz1-tmp/ST_DIVERSITY/PUB/USER/zhouyang/bin/convert/personal/';
use GFF;

die "perl $0 <gff> <genome.fa>" unless @ARGV == 2;
my $fasta = &read_fasta($ARGV[1]);

my ($in_file, $ref) = @_;
my %tmp;
open(IN, $ARGV[0]) or die $!;
while (<IN>) {
	my @info = split /\t/;
	next unless $info[2] eq "CDS";
	die unless $info[8] =~ /^Parent=(\S+?);/;
	my $id = $1;
	push @{$tmp{$id}}, [$info[0], $info[3], $info[4], $info[6]]; ## chr bg ed strand
}
close IN;

foreach my $id (sort keys %tmp) {
	my $strand = $tmp{$id}->[0]->[3];
	if ($strand eq "+") {
		@{$tmp{$id}} = sort {$a->[1] <=> $b->[1]} @{$tmp{$id}};
	} elsif ($strand eq "-") {
		@{$tmp{$id}} = sort {$b->[1] <=> $a->[1]} @{$tmp{$id}};
	} else {
		die;
	}
	my ($start, $end);
	my @range;
	for(my $i=0;$i<@{$tmp{$id}};$i++){
		my ($chr, $bg, $ed, $strand) = @{$tmp{$id}[$i]};
		if($strand eq '+'){
			for(my $j=$bg;$j<=$ed;$j++){
				push @range, $j;
			}
		}
		else{
			for(my $j=$ed;$j>=$bg;$j--){
				push @range, $j;
			}
		}
	}
	my $chr = $tmp{$id}[0][0];
	next  if (@range % 3 != 0);
	for(my $i=0;$i<@range;$i+=3){
		my $idx = $i/3+1;
		my $b1 = substr($$fasta{$chr}, $range[$i]-1, 1);
		my $b2 = substr($$fasta{$chr}, $range[$i+1]-1, 1);
		my $b3 = substr($$fasta{$chr}, $range[$i+2]-1, 1);
		if ($range[$i+1] < $range[$i]){
			# strand = '-'
			$b1 = complement($b1);
			$b2 = complement($b2);
			$b3 = complement($b3);
		}
		my $codon = $b1.$b2.$b3;
		my $amino = &translate($codon);
		print "$id\t$chr\t$strand\t$range[$i]:$b1\t$range[$i+1]:$b2\t$range[$i+2]:$b3\t$codon\t$amino\t$idx\n";
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

sub read_fasta{
	my $file = shift;
	my %fasta;
	open IN, $file or die "$!";
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		my $head = $_;
		chomp $head;
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		$seq =~ s/\n//g;
		$/="\n";
		$fasta{$head} = $seq;
	}
	return \%fasta;
}

sub complement{
	my $base = shift;
	$base=~tr/AGCTagct/TCGAtcga/;
	return $base;
}



