#!/usr/bin/perl

=head1  Name

  to.link.pl

=head1  Vesion

  Authot:yinxinxin  yinxinxin@genomics.cn

=head1  Description

  This procedure is used to connect two file to one according to the rows .
  You can chose no blank bewteen two rows .defalut no blank.
  Note that the lines of two files are likely different .So we hope you can use "--strict yes" all the time. 
=head1  Usage
      
  perl to.link.pl <file1> <file2> ... [options]
  (The files are sequential.)

         --blank     yes || no .To chose whether there are blank between two linked row or not.default "no".
         --help      show  help information.     
         --strict    yes || no .To set whether allowing different  line numbers of two files or not.default "no".
      
=head1  Example

  perl to.link.pl file1 file2 
  perl to.link.pl file1 file2 --strict yes --blank yes 
  
  file1:                              file2:                                linked file:  perl to.link.pl file1 file2 
        afsadfsadfsadfsadf                  FADSFDSADFSDAFSADFSA                      afsadfsadfsadfsadfFADSFDSADFSDAFSADFSA
        afsadfsadfsadfsfsdafsa              AFSADFSADFSDFSADFSADF                     afsadfsadfsadfsfsdafsaAFSADFSADFSDFSADFSADF


=cut

use strict;
use warnings;
use Getopt::Long;
my ($Blank,$Help,$Strict);

GetOptions(
          "blank:s"=>\$Blank,
	  "help!"=>\$Help,
	  "strict:s"=>\$Strict,
	  )
;
die `pod2text $0` if(@ARGV==0 || $Help);
$Blank ||= "no";
$Strict ||="no";
my($file1,$file2);
open IN1,shift;
my(@t1,@t2,@tmp_all);
@t1=<IN1>;
close IN1;
my $decide=1;
while($decide){
	  $file2=shift;
          last if ( !defined $file2) ;
	  open IN2,$file2;
	  my @t2=<IN2>;
	  my ($max,$min);my $space;
	  my $number1=scalar(@t1);my $number2=scalar(@t2);
	  if ( $Strict=~/[Yy][Ee][Ss]/){
	  die "Warnings : the lines of two files are different" if ( $number1!=$number2);
	  }

	  if ($number1>$number2){
               $max=$number1;
	       $min=$number2;
	       $space=""; 
	    }
	  else{
	       $max=$number2;
	       $min=$number1;
	       my $lengthpre=length($t1[$min-1]);
	       $space="\t"." "x$lengthpre;
	       if($Blank=~/[Yy][Ee][Ss]/){ $space="\t"." "x$lengthpre;}
	       else{$space=""." "x$lengthpre;}
	    }

	  if ($Blank=~/[nN][Oo]/){
	      for (my $i=0;$i<$max;$i++){
	         if($i<$min){
	         chomp($t1[$i]);chomp($t2[$i]);
                 $tmp_all[$i]=$t1[$i].$t2[$i];
	         }
	         else{
	           if (defined $t1[$i]){chomp($t1[$i]);$tmp_all[$i]=$space.$t1[$i];}
	           if (defined $t2[$i]){chomp($t2[$i]);$tmp_all[$i]=$space.$t2[$i];}
	         }
	      }
	  }
	  elsif($Blank=~/[Yy][Ee][Ss]/){
	       for (my $i=0;$i<$max;$i++){
	           if($i<$min){
	           chomp($t1[$i]); chomp($t2[$i]);
	           $tmp_all[$i]=$t1[$i]."\t".$t2[$i];
	           }
	           else{
	             if (defined $t1[$i]){chomp($t1[$i]);$tmp_all[$i]=$space.$t1[$i]; }
	             if (defined $t2[$i]){chomp($t2[$i]);$tmp_all[$i]=$space.$t2[$i]; }
		     }
               }
	  }
	  @t1=@tmp_all;
	  close IN2;
          } 

foreach(@tmp_all){
  print $_,"\n";
  }
