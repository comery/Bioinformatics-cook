#!/usr/bin/perl -w
# written by BGI-tech at April 2013
use strict;
use Getopt::Long;
my ($input,$output,@tmp,$chr,$pos,$ref,$alt,@allele,$help);

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
	perl $0 -input pear.fina.vcf -output pear.final.snp.genotype
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
        my @tmp = split ;
        my $chr = $tmp[0];
        my $pos = $tmp[1];
        my $ref = $tmp[3];
        my $alt = $tmp[4];
                       
        my $output = $chr."\t".$pos."\t".$ref."\t";
        for (my $sample = 9;$sample <=$#tmp ;$sample ++ ){
                my @arr = split ':' , $tmp[$sample]; 
                if ($alt=~/^\w$/) {
                if ($tmp[$sample] eq "./." ){
                        $output.="-";
	
                }
                elsif ($arr[0] eq "0/0"){
                        $output.=$ref;
                }
                elsif ($arr[0] eq "1/1"){
                        $output .=$alt;
                }
                else{
                        $output .=$combin{$ref.$alt};
                }
              }            
              if ($alt=~/,/) {  
              	@allele=split /,/,$alt;            	
              	if ($tmp[$sample] eq "./." ){
                        $output.="-";
                }
                elsif ($arr[0] eq "0/0"){
                        $output.=$ref;
                }
                elsif ($arr[0] eq "1/1"){
                        $output .=$combin{$allele[0].$allele[0]};
                }
                 elsif ($arr[0] eq "2/2"){
                        $output .=$combin{$allele[1].$allele[1]};
                        }
                 elsif ($arr[0] eq "0/1"){
                        $output .=$combin{$ref.$allele[0]};
                        }
                elsif ($arr[0] eq "0/2"){
                        $output .=$combin{$ref.$allele[1]};
                }
                elsif ($arr[0] eq "1/2"){
                        $output .=$combin{$allele[0].$allele[1]};
                }                    
                
              }         
                if ($sample < $#tmp){
                        $output .=" ";
                }
        }
        print OUTPUT $output."\n";
}


