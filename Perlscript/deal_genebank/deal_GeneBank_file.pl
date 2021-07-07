#! /usr/bin/perl -w
=head1 Description

	this script is used to find 18s/ITS1/5.8S/ITS2 and 28S ribosomal RNA information from genebank file,and extract GI number and corresponding sequence ; 

=head1 Usage
		
	perl deal_GeneBank_file.pl <.gb>

=head1 Version
		
	Version :  1.0
	Created :  06/25/2014 09:28:12 PM

=head1 Contact

	Author       :	Shiyi Du
	Email        :	dushiyi@genomics.cn
	Organization :	BGI-NGB 

=cut

use strict;
use warnings;
#my $file1=shift;     #genebank file
my ($version,@GI,$gi,$panduanRNA,$panduanITS,$panduan,$origin,$len,$fa,$panduan1,$panduan2,$panduan3,$panduan4,$panduan5,$panduan6,$panduan7);
my ($F1,$E1,$C1,$F2,$E2,$C2,$F3,$E3,$C3,$F4,$E4,$C4,$F5,$E5,$C5,$F6,$F7,$E6,$E7,$C6,$C7,$S);

#open GB,"<$file1" or die;
open OUT1,">NCBI_VCP_genus_GeneBank_18S_rRNA.fasta";
open OUT2,">NCBI_VCP_genus_GeneBank_5.8S_rRNA.fasta";
open OUT3,">NCBI_VCP_genus_GeneBank_28S_rRNA.fasta";
open OUT4,">NCBI_VCP_genus_GeneBank_ITS1_rRNA.fasta";
open OUT5,">NCBI_VCP_genus_GeneBank_ITS2_rRNA.fasta";
die `pod2text $0` unless (@ARGV == 1);

while($S||=<>){
         chomp $S;
         my $line=$S;
	 $S=();
	 $line=~s/^\s+//g;
         if($line=~/^VERSION(.+)/){
                  $version=$1;
                  @GI=split /:/,$version;
                  $gi=$GI[1];
                  $gi=~s/\s//g;
		  next;
         }	 
	 if($line=~/rRNA\s+\<?(\d+)..\>?(\d+)/){
		  my $first_number=$1;
                  my $end_number=$2;
		  while ($panduanRNA=<>){
			   if ($panduanRNA=~/\//){
				    if ($panduanRNA=~/\/.+="18S/){
					     $panduanRNA=~s/\s//g;
					     $panduan1=$panduanRNA;
					     $F1=$first_number;
					     $E1=$end_number;
					     $C1=$E1-$F1+1;
					     last;
				    }
				    if ($panduanRNA=~/\/.+="5.8S/){
					     $panduanRNA=~s/\s//g;
					     $panduan2=$panduanRNA;
					     $F2=$first_number;
					     $E2=$end_number;
					     $C2=$E2-$F2+1;
					     last;
				    }
				    if ($panduanRNA=~/\/.+="28S/) {
					     $panduanRNA=~s/\s//g;
					     $panduan3=$panduanRNA;
					     $F3=$first_number;
					     $E3=$end_number;
			   		     $C3=$E3-$F3+1;
					     last;
				    }
			   	    next;
			   }
			   else {
				    $S=$panduanRNA;
				    last;
			   }
			   next;
		  }              
		  next;
         }
         if ($line=~/misc_RNA\s+\<?(\d+)..\>?(\d+)/) {
                  my $first_number=$1;
                  my $end_number=$2;
                  my $c=$end_number-$first_number;
                  $panduanITS=<>;
                  if ($panduanITS=~/\/.+="internal\stranscribed\sspacer\s1/) {
                           $panduanITS=~s/\s//g;
                           $panduan4=$panduanITS;
			   $F4=$first_number;
			   $E4=$end_number;
			   $C4=$E4-$F4+1;
			   next;
                  }
                  if ($panduanITS=~/\/.+="internal\stranscribed\sspacer\s2"/) {
                           $panduanITS=~s/\s//g;
			   $panduan5=$panduanITS;
			   $F5=$first_number;
			   $E5=$end_number;
			   $C5=$E5-$F5+1;
			   next;
                  }
		  next;
         }
	 
	 if ($line=~/misc_feature\s+\<?(\d+)..\>?(\d+)/) {
		  my $first_number=$1;
                  my $end_number=$2;
                  my $c=$end_number-$first_number;
                  $panduan=<>;
		  if ($panduan=~/\/note="internal\stranscribed\sspacer\s1/) {
                           $panduan=~s/\s//g;
			   $panduan6=$panduan;
                           $F6=$first_number;
			   $E6=$end_number;
			   $C6=$E6-$F6+1;
			   next;
		  }
		  if ($panduan=~/\/note="internal\stranscribed\sspacer\s2/) {
                           $panduan=~s/\s//g;
                           $panduan7=$panduan;
			   $F7=$first_number;
			   $E7=$end_number;
			   $C7=$E7-$F7+1;
			   next;
                  }
		  next;
	 }
	 if ($line=~/ORIGIN/) {
		  $/="//";
                  my $list=<>;
                  chomp $list;
		  $/="\n";
                  $list=~ s/\s+//g;
                  $list=~ s/\d//g;
                  $list=~ tr/[a-z]/[A-Z]/;
                  $len=length ($list);
		  if (defined $panduan1) {
			   if ($C1>100) {
				    my $fa1= substr ($list,$F1,$C1);
				    print OUT1 ">gi:$gi:$F1:$E1:$C1:$panduan1 $len\n$fa1\n";
			   }else {
				    next;
			   }
		  }
		  if (defined $panduan2) {
			   if ($C2>100) {
				    my $fa2= substr ($list,$F2,$C2);
				    print OUT2 ">gi:$gi:$F2:$E2:$C2:$panduan2 $len\n$fa2\n";
			   }else{
				    next;
			   }
		  }
		  if (defined $panduan3) {
			   if ($C3>100) {
				    my $fa3= substr ($list,$F3,$C3);
				    print OUT3 ">gi:$gi:$F3:$E3:$C3:$panduan3:$len\n$fa3\n";
			   }else{
				    next;
			   }
		  }
		  if (defined $panduan4) {
			   if ($C4>100) {
				    my $fa4= substr ($list,$F4,$C4);
				    print OUT4 ">gi:$gi:$F4:$E4:$C4:$panduan4:$len\n$fa4\n";
			   }else{
				    next;
			   }
		  }
		  if (defined $panduan5) {
			   if ($C5>100) {
				    my $fa5= substr ($list,$F5,$C5);
				    print OUT5 ">gi:$gi:$F5:$E5:$C5:$panduan5:$len\n$fa5\n";
			   }else{
				    next;
			   }
		  }
		  if (defined $panduan6) {
			   if ($C6>100) {
				    my $fa6= substr ($list,$F6,$C6);
				    print OUT4 ">gi:$gi:$F6:$E6:$C6:$panduan6:$len\n$fa6\n";
			   }else{
				    next;
			   }
		  }
		  if (defined $panduan7) {
			   if ($C7>100) {
				    my $fa7= substr ($list,$F7,$C7);
				    print OUT5 ">gi:$gi:$F7:$E7:$C7:$panduan7:$len\n$fa7\n";
			   }else{
				    next;
			   }
		  }
	 }else{
		  next;
	 }
         next;
}
#close GB;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;