#!/usr/bin/perl -w 

=head1 Name
find the Mutation frequency of synonymous and non-synonymous mutations for each site

=head1 Description
 The output file format is like this :
 	site	deepth	consensus	befor_mut	after_mut	site_info	fruency
 	28      3319    A       ATG => M        TTG => L        A => T  0.000301295570955107

=head1 Contact & Version
  Author: Chentao YANG, yangchentao@genomics.org.cn
  Version: 1.0,  Date: 2014-8-8

=head1 Command-line Option
  perl   <infile | STDIN>
  --start <number>           the number of CDS start site
  --end   <number>           the number of CDS	end   site
  --pro   <string>           the protein this domain coded
  --help                     output help information to screen  

=head1 Usage Exmples
  perl  ./syn_and_nsyn.pl  Influenza-A-seg1 -start  1 -end 2280  
  
=cut

#################################################################################
#use strict;
use Getopt::Long;
use List::Util qw(max);

my ($start,$end,$Help,$cds_len,$inputfile,$outputdir,$protein);
GetOptions(
	"start:n"=>\$start,
	"end:n"=>\$end,
	"help"=>\$Help,
	"pro:s"=>\$protein,
	);
die `pod2text $0` if (@ARGV==0 || $Help);
#`mkdir $outputdir` if (! -e "./$outputdir");
$inputfile = $ARGV[0];
my $segment = $1 if ($inputfile =~ /outdir\/(seg\d)\/\w+/); ## outdir name is based on inputfile name 
#print "$segment\n";
$outputdir = "outdir/"."$segment";
$cds_len=$start-$end;
open IN,"<$inputfile" or die "please give me a right file pathway!";
#print "$start,$end";

############################ localize the start place and make a hash for querry
open MUT,">$outputdir/mut_rate.txt";
my (%sites,@sub_sites,@lines,%deepth,%count,$line,$len,$key,%key_base,$start_key,$min_key);
while (<IN>){
	my @info=split /\s+/,$_;
	$key=$info[1]; ## this key is site number in inputfle
	my $aa="A";my $t="T";my $g="G";my $cc="C";
	push @sub_sites,$aa if ($info[3]>0);
	push @sub_sites,$t if ($info[4]>0);
	push @sub_sites,$g if ($info[5]>0);
	push @sub_sites,$cc if ($info[6]>0);
	$deepth{$key}=$info[2];
	$sites{$key}="@sub_sites";
	$count{$key}{$aa}=$info[3];
	$count{$key}{$t}=$info[4];
	$count{$key}{$g}=$info[5];
	$count{$key}{$cc}=$info[6];
	$len=$key;
	@sub_sites=();
	my %base_num=(
					$info[3] => $aa,
					$info[4] => $t,
					$info[5] => $g,
					$info[6] => $cc,
					);
	my $max_num=max(keys %base_num); ## the deepest site is orginal site
	$key_base{$key}=$base_num{$max_num};

	my $mut_rate = 1-($max_num/$info[2]);
	print MUT "$key\t$mut_rate\n";
}
	my @sorted_keys=sort{$a<=>$b}keys %key_base;
	$min_key=$sorted_keys[0];
	my @useless=splice(@sorted_keys,0,$start-$min_key);
	#print "querry__\t";
foreach (@sorted_keys){
    #print "$key_base{$_}";
}

## standard codon table
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

open OUT1,">$outputdir/$protein\_nsyn.xls";
open OUT2,">$outputdir/$protein\_syn.xls";
my ($j,$i);
for ($i=$start;$i<=($end-2);$i+=3) {
	$j=($i-$start)/3+1;
	my $a=$key_base{$i};
	my $b=$key_base{$i+1};
	my $c=$key_base{$i+2};
	my $ref_thr="$a"."$b"."$c";
	my $ref_pro=$code{$ref_thr};
	#print "$ref_thr";
	foreach my $x(split /\s+/,$sites{$i}) {
		my $tem_thr1="$x"."$b"."$c";
		if ($code{$tem_thr1} ne $ref_pro ) {
			my $syn1=&rate($count{$i}{$x},$deepth{$i});
			print OUT1 "$protein\t$j\t$i\t$deepth{$i}\t$key_base{$i}\t$ref_thr => $ref_pro\t$tem_thr1 => $code{$tem_thr1}\t$a => $x\t$syn1\n";
		}elsif($code{$tem_thr1} eq $ref_pro && $ref_thr ne $tem_thr1 ) {
			my $nsyn1=&rate($count{$i}{$x},$deepth{$i});
			print OUT2 "$protein\t$j\t$i\t$deepth{$i}\t$key_base{$i}\t$ref_thr => $ref_pro\t$tem_thr1 => $code{$tem_thr1}\t$a => $x\t$nsyn1\n";
		}
	}
	foreach my $y(split /\s+/,$sites{$i+1}) {
		my $tem_thr2="$a"."$y"."$c";
		my $ii=$i+1;
		if ($code{$tem_thr2} ne $ref_pro ) {
			my $syn2=&rate($count{$ii}{$y},$deepth{$ii});
			print OUT1 "$protein\t$j\t$ii\t$deepth{$ii}\t$key_base{$ii}\t$ref_thr => $ref_pro\t$tem_thr2 => $code{$tem_thr2}\t$b => $y\t$syn2\n";
		}elsif($code{$tem_thr2} eq $ref_pro && $ref_thr ne $tem_thr2 ) {
			my $nsyn2=&rate($count{$ii}{$y},$deepth{$ii});
			print OUT2 "$protein\t$j\t$ii\t$deepth{$ii}\t$key_base{$ii}\t$ref_thr => $ref_pro\t$tem_thr2 => $code{$tem_thr2}\t$b => $y\t$nsyn2\n";
		}
	}
	foreach my $z(split /\s+/,$sites{$i+2}) {
		my $tem_thr3="$a"."$b"."$z";
		my $iii=$i+2;
		if ($code{$tem_thr3} ne $ref_pro ) {
			my $syn3=&rate($count{$iii}{$z},$deepth{$iii});
			print OUT1 "$protein\t$j\t$iii\t$deepth{$iii}\t$key_base{$iii}\t$ref_thr => $ref_pro\t$tem_thr3 => $code{$tem_thr3}\t$c => $z\t$syn3\n";
		}elsif($code{$tem_thr3} eq $ref_pro && $ref_thr ne $tem_thr3 ) {
			my $nsyn3=&rate($count{$iii}{$z},$deepth{$iii});		
			print OUT2 "$protein\t$j\t$iii\t$deepth{$iii}\t$key_base{$iii}\t$ref_thr => $ref_pro\t$tem_thr3 => $code{$tem_thr3}\t$c => $z\t$nsyn3\n";
		}
		}
	}
close IN;
close OUT1;
close OUT2;
		
			
##############sub#######################
sub rate {
	my @array=@_;
	my $rate=$array[0]/$array[1];
	return $rate;
	}
########################################







