#!/usr/bin/perl
=head1 Name
statistic the ENC(effective number of codon usage) with a new method (Sun X, Yang Q, Xia X. An improved implementation of effective number of codons (Nc)[J]. Molecular biology and evolution, 2013, 30(1): 191-196.)

=head1 Description
The effective number of codons (Nc) is a widely used index for characterizing codon usage bias because it does not require a set of reference genes as does codon adaptation index (CAI) and because of the freely available computational tools such as CodonW. However,Nc, as originally formulated has many problems. For example, it can have values far greater than the number of sense codons; it treats a 6-fold compound codon family as a single-codon family although it is made of a 2-fold and a 4-fold codon family that can be under dramatically different selection for codon usage bias; the existing implementations do not handle all different genetic codes; it is often biased by codon families with a small number of codons. We developed a newNcthat has a number of advantages over the original Nc. Its maximum value equals the number of sense codons when all synonymous codons are used equally, and its minimum value equals the number of codon families when exactly one codon is used in each synonymous codon family. It handles all known genetic codes. It breaks the compound codon families (e.g., those involving amino acids coded by six synonymous codons) into 2-fold and 4-fold codon families. It reduces the effect of codon families with few codons by introducing pseudocount and weighted averages. The new Nc has significantly improved correlation with CAI than the original Nc from CodonW based on protein-coding genes from Saccharomyces cerevisiae, Caenorhabditis elegans,Drosophila melanogaster, Escherichia coli, Bacillus subtilis, Micrococcus luteus,andMycoplasma genitalium.It also correlates better with protein abundance data from the yeast than the original Nc.

=head1 Contact & Version
  Author: Chentao YANG, yangchentao@genomics.org.cn
  Version: 1.0,  Date: 2014-11-23

=head1 Command-line Option
  perl   <infile | STDIN>
  --in    <string>           the seq of gene           
  --out   <string>           default STDOUT
  --help                     output help information to screen  

=head1 Usage Exmples
  perl  new_ENC.pl  -in  seq.fa -out dir

=head1 Output file format like this:
###T2D-6A_GL0083352  [gene]  locus=scaffold27901_4:37993:126222:-[Complete]     29410
F       TTT=>346        0.0117647058823529      TTC=>248        0.00843250595035702     594     0.0540049095372307      6.88964976516282e-07	3.72074912310821e-08
note:	first line is gene id and the number of amino;
		second line is "amino	codon1	frequency_codon1	... 	sum_number_all_codon	precentage of i-fold all codons		Fcf value"
  
=cut
#################################################################################
use strict;
use warnings;
use Getopt::Long;
my ($Help,$inputfile,$outdir);
GetOptions(
	"help" => \$Help,
	"out:s" => \$outdir
	);
$outdir||= "./";
die `pod2text $0` if (@ARGV==0);
die `pod2text $0` if ($Help);
$inputfile = shift;
open IN,"$inputfile" or die "please give me a right file pathway!";
open OUT,">$outdir/$inputfile.ENC.xls";
############################ Code Table #############################################

my %code =(
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
			);
########################### Main ######################################################
my ($id,$seq,$length,$amino_len,%codon_usage,$key,$val,$num_syn);

$/=">";<IN>;$/="\n";
while (<IN>){
	$id = $_;
	chomp $id;
	$/=">";
	$seq = <IN>;
	chomp $seq;
	$seq =~ s/\n//g;
	$length = length($seq);
	$amino_len=$length/3;
	while (($key,$val) = each %code) {
		$codon_usage{$val}{$key} = 0;# estiblish a hash table
	}
	for (my $i = 0;$i <= $length-1;$i+=3) {
		my $tmpstr = substr($seq,$i,3);
		if (defined $codon_usage{$code{$tmpstr}}{$tmpstr}) {
			$codon_usage{$code{$tmpstr}}{$tmpstr}++;
		}else{
			$codon_usage{$code{$tmpstr}}{$tmpstr} = 1;
		}
	}
########################## cycle first time #############################################
	my $num_k1 =0;
	my $num_k2 =0;
	my $num_k3 = 0;
	my $num_k4 = 0;
	my $num_k6 = 0;
	my $TER = 0;
	my %KI_count;
	print OUT "###$id\t$amino_len\n";
	foreach my $key1 (keys %codon_usage) {
		my $hash2 = $codon_usage{$key1};
		$num_syn = scalar (keys %$hash2);			#the number of synonymous codons for a amino
		foreach my $key2 (sort{$hash2->{$b}<=>$hash2->{$a}} keys %$hash2){ #scan each synonmous codon 
			if ($num_syn == 1) {$num_k1 +=  $hash2->{$key2};}
			if ($num_syn == 2) {$num_k2 +=  $hash2->{$key2};}
			if ($num_syn == 3 && $key1 ne "U") {$num_k3 +=  $hash2->{$key2};}
			if ($num_syn == 4) {$num_k4 +=  $hash2->{$key2};}
			if ($num_syn == 6) {$num_k6 +=  $hash2->{$key2};}
			if ($key1 eq "U")  {$TER += $hash2->{$key2};}
		}
		print STDERR  "there are $TER terminal codons,Be careful!" if ($TER >1);
	}
	
		%KI_count =(		1 => "$num_k1",
							2 => "$num_k2",
							3 => "$num_k3",
							4 => "$num_k4",
							6 => "$num_k6"
					);

######################### cycle second time ###############################################
	my ($FCF_amino,$nj_multi_Fcf,@K2_nF,@K3_nF,@K4_nF,$NJ);
	foreach my $key1 (keys %codon_usage) {
		my $hash2 = $codon_usage{$key1};
		$num_syn = scalar (keys %$hash2);			#the number of synonymous codons for a amino
		print OUT "$key1\t";						#print this amino
		my $sum_num_a_amino = 0;
		foreach my $key2 (sort{$hash2->{$b}<=>$hash2->{$a}} keys %$hash2){ #scan each synonmous codon 
			my $frequency = &Ni($hash2->{$key2});	#one codon's frequency of appearance
			$FCF_amino = &Fcf(values %$hash2);		#the value of Fcf(see the detail in mentioned refrence)
			print OUT $key2."=>".$hash2->{$key2}."\t".$frequency."\t";
			$sum_num_a_amino += $hash2->{$key2};	#sum of each codon's frequency in a codon family
		}
		$NJ = $sum_num_a_amino/$KI_count{$num_syn};	#nj
		$nj_multi_Fcf = $NJ*$FCF_amino;				#nj*Fcf
		print OUT "$sum_num_a_amino\t$NJ\t$FCF_amino\t$nj_multi_Fcf\n";
		if ($num_syn == 2) {push @K2_nF,$nj_multi_Fcf;}
		if ($num_syn == 3 && $key1 ne "U") {push @K3_nF,$nj_multi_Fcf;}
		if ($num_syn == 4) {push @K4_nF,$nj_multi_Fcf;}
	}
#	print "@K3_nF\n";
	my ($k2_nf,$k3_nf,$k4_nf,$Nc);
	$k2_nf = &sum(@K2_nF);
	$k3_nf = &sum(@K3_nF);
	$k4_nf = &sum(@K4_nF);
	my ($k2,$k3,$k4);
	$k2 = 9/$k2_nf;
	$k3 = 1/$k3_nf;
	$k4 = 5/$k4_nf;
	$Nc = 2+$k2+$k3+$k4;
	print "2+$k2+$k3+$k4 = $Nc\n";

	$/="\n";
}
close IN;
close OUT;

########################### sub scripts ###################################
sub Ni {
	my $a = shift;
	my $ni = $a/$amino_len;
	return $ni
}
##############################
sub Fcf {
	my @array = @_;
	my $m = $num_syn;
	my $n = $amino_len;
	my $fcf;
	foreach my $i(@array) {
		my $squ = ($i+1)/($m+$n)**2;# here you can modif "my $squ = (&Ni($i)+1)/($m+$n)**2";that means  codon number rate in all codons
		$fcf += $squ;
	}
	return $fcf;
}
#############################
sub sum {
	my @array1 = @_;
	my $sum = 0;
	foreach my $q(@array1){
		$sum += $q;
	}
return $sum;
}
############################
