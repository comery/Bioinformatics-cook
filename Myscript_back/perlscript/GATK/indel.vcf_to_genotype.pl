#!/usr/bin/perl 
# written by BGI-tech at April 2013
use strict;
use Getopt::Long;

my ($input,$output,@tmp,$chr,$pos,$ref,$alt,%tmp,@allele,$help);

GetOptions("input:s" => \$input, "output:s" => \$output, "help|?" => \$help);

if (!defined $input ||  defined $help) {
	die << "USAGE";
description: from vcf file get genotypefile
usage: perl $0 [options]
options:
	-input <str> *  input vcf file ,may sigle or multiple sample;
	-output <str>*  output genotype file ;
  -help|?         help information;
e.g.:
	perl $0 -input pear.fina.vcf -output pear.final.indel.genotype
USAGE
}

my %combin = (
		"AC"=>   "M" , 
		"AG"=>   "R" , 
		"AT"=>   "W" , 
		"CT"=>   "Y" , 
		"CG"=>   "S" , 
		"GT"=>   "K" , 

		"AA"=>   "A" , 
		"TT"=>   "T" , 
		"GG"=>   "G" , 
		"CC"=>   "C" , 

		"CA"=>   "M" , 
		"GA"=>   "R" , 
		"TA"=>   "W" , 
		"TC"=>   "Y" , 
		"GC"=>   "S" , 
		"TG"=>   "K" , 
	     );
open VCF , $input or die $!;
open OUTPUT , ">$output" or die $!;

while (<VCF>){
	next if /^#/;
	my @tmp = split /\t/,$_ ;
	my $chr = $tmp[0];
	my $pos = $tmp[1];
	my $ref = $tmp[3];
	my $alt = $tmp[4];
	
	my $output = $chr."\t".$pos."\t".$ref."\t";
	for (my $sample = 9;$sample <=$#tmp ;$sample ++ ){
		my  @arr = split ':' , $tmp[$sample];
		if ($tmp[$sample] eq "./." ){
                        $output.="-";
}
		elsif (($arr[0] eq "0/1" &&  $arr[2] < 5) ||$arr[2] < 3){
			$output .="-";
			}
		elsif ($arr[0] eq "0/1") {
		$output .="$tmp[3]/$tmp[4]";
	}
	elsif ($arr[0] eq "0/0") {
		$output .="$tmp[3]/$tmp[3]";
	}
	elsif ($arr[0] eq "1/1") {
		$output .="$tmp[4]/$tmp[4]";
	}
		
		
		if ($sample < $#tmp){
			$output .=" ";
		}
	}
	
	print OUTPUT  $output."\n";
}


