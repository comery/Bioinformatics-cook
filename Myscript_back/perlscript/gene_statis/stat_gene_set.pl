#!/usr/bin/perl
use strict;
use warnings;

die "Usage:$0 <*gff|*gff2> \n" if @ARGV<1;

my $gff=shift;

my %Gene;
open IN,$gff or die "$!";
while(<IN>){
	next if /^#/;
	chomp;
	my @c=split(/\t/);
	@c[3,4]=@c[4,3] if $c[3]>$c[4]; 
	if ($c[2] eq 'mRNA' && ( $c[8]=~/ID=([^;]+)/ || $c[8]=~/GenePrediction\s+(\S+)/)){
		@{$Gene{$1}{mRNA}}=@c;
	}elsif($c[2] eq 'CDS' && ( $c[8]=~/Parent=([^;]+)/ || $c[8]=~/GenePrediction\s+(\S+)/)){
		push @{$Gene{$1}{CDS}},[@c];
	}
}
close IN;

my %Stat;

foreach my $id ( sort keys %Gene ){
	$Stat{gene_num}++;
	$Stat{gene_len}+=abs($Gene{$id}{mRNA}[4]-$Gene{$id}{mRNA}[3])+1;
	@{$Gene{$id}{CDS}}=sort {$a->[3]<=>$b->[3]} @{$Gene{$id}{CDS}};
	for(my $i=0;$i<@{$Gene{$id}{CDS}};$i++){
		$Stat{exon_num}++;
		$Stat{exon_len}+=abs($Gene{$id}{CDS}[$i][4]-$Gene{$id}{CDS}[$i][3])+1;
		next if $i==0;
		$Stat{intron_num}++;
		$Stat{intron_len}+=$Gene{$id}{CDS}[$i][3]-$Gene{$id}{CDS}[$i-1][4]-1;
	}
}

$Stat{ave_gene}=$Stat{gene_len}/$Stat{gene_num};
$Stat{ave_cds}=$Stat{exon_len}/$Stat{gene_num};
$Stat{ave_exon_num}=$Stat{exon_num}/$Stat{gene_num};
$Stat{ave_exon}=$Stat{exon_len}/$Stat{exon_num};
$Stat{ave_intron}=$Stat{intron_len}/$Stat{intron_num};

print "gene_number\taverage_gene_len\taverage_cds_len\taverage_exon_number\taverage_exon_len\taverage_intron_len\n";
printf "%.2f\t",$Stat{gene_num};
printf "%.2f\t",$Stat{ave_gene};
printf "%.2f\t",$Stat{ave_cds};
printf "%.2f\t",$Stat{ave_exon_num};
printf "%.2f\t",$Stat{ave_exon};
printf "%.2f\t",$Stat{ave_intron};
print "\n";
