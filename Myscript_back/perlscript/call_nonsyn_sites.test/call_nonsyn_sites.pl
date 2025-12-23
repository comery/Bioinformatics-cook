#!usr/bin/perl -w 
#use strict;
use List::Util qw(max);
open IN,shift;
my (%sites,@sub_sites,@lines,$line,$len,$key,%key_base,$start_key,$end_key);
############################ localize the start place and make a hash for querry
while (<IN>){
	my @info=split /\s+/,$_;
	$key=$info[1];
	my $aa="A";my $t="T";my $g="G";my $cc="C";
	push @sub_sites,$aa if ($info[3]>0);
	push @sub_sites,$t if ($info[4]>0);
	push @sub_sites,$g if ($info[5]>0);
	push @sub_sites,$cc if ($info[6]>0);
	$sites{$key}="@sub_sites";
	$len=$key;
	@sub_sites=();
	my %base_num=(
					$info[3] => $aa,
					$info[4] => $t,
					$info[5] => $g,
					$info[6] => $cc,
					);
	my $max_num=max(keys %base_num);
	$key_base{$key}=$base_num{$max_num};
}
my @sorted_keys=sort{$a<=>$b}keys %key_base;
foreach $k(@sorted_keys){
#print "$k\n";
	my $thr_base="$key_base{$k}"."$key_base{$k+1}"."$key_base{$k+2}";
	if ($thr_base eq "ATG"){
	$start_key=$k;
	last;
	}
}
#for($k=$start_key;$k<=($len/3);$k++){
#	my  $thr_base1="$key_base{$k}"."$key_base{$k+1}"."$key_base{$k+2}";
#	if  ($thr_base1 eq "TAA"||$thr_base1 eq "TGA"||$thr_base1 eq "TAG") {
#		$end_key=$k;
#		print "$end_key\n";
#		last;
#	}
	
#}
############################ extract the complete cds and make a hash for refrence
open IN1, shift;
my ($ref_line,%ref_sites,$ref_line_old,$cds_len);
$/=">";<IN1>;$/="\n";
my $id=<IN1>;
my @ids=split /\s+/,$id;
my $start=$ids[3];
my $end=$ids[4];
$cds_len=$end-$start+1;

while (<IN1>){
	chomp;
	$ref_line_old.=$_;
	$ref_line=substr($ref_line_old,$start-1,$cds_len);
	}
	
for (my $j=1;$j<=$cds_len;$j++) {
	my $ref_site=substr $ref_line,$j-1,1;
	$ref_sites{$j}=$ref_site;
	}


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

for ($i=$start_key,$y=1;$i<=($len/3),$y<=($cds_len/3);$i+=3,$y+=3) {
	my $a=$ref_sites{$y};
	my $b=$ref_sites{$y+1};
	my $c=$ref_sites{$y+2};
	my $ref_thr="$a"."$b"."$c";
	my $ref_pro=$code{$ref_thr};
	foreach my $d(split /\s+/,$sites{$i}) {
		my $tem_thr1="$d"."$b"."$c";
		print "$i\t$ref_thr => $ref_pro\t$tem_thr1 => $code{$tem_thr1}\t$a => $d\n" if ($code{$tem_thr1} ne $ref_pro ||$ref_thr ne $tem_thr1 );
		}
	foreach $d(split /\s+/,$sites{$i+1}) {
		my $tem_thr2="$a"."$d"."$c";
		my $ii=$i+1;
		print "$ii\t$ref_thr => $ref_pro\t$tem_thr2 => $code{$tem_thr2}\t$b => $d\n" if ( $code{$tem_thr2} ne $ref_pro ||$ref_thr ne $tem_thr2);
		}
	foreach $d(split /\s+/,$sites{$i+2}) {
		my $tem_thr3="$a"."$b"."$d";
		my $iii=$i+2;
		print "$iii\t$ref_thr => $ref_pro\t$tem_thr3 => $code{$tem_thr3}\t$c => $d\n" if ( $code{$tem_thr3} ne $ref_pro || $ref_thr ne $tem_thr3);
		}
	}
close IN;
close IN1;
		
			
















