#! /usr/bin/perl

use strict;
use warnings;

=head1
Description

     This procedure is used to simplify the process of qsub.
     you can modify it according to you needs.

Usage 
     
     perl help.qsub.pl  [option] <commonder>
 
                     --help           show the help information.
		     --l              to set the memory of the procedure.
                     --shell          yes|no, to confirm the commander is  shell or not. if there is no input string about shell, then the procedure will check automatically
                     
=cut


my ($Help,$Shell,$L);
use Getopt::Long;
GetOptions(
	   "help!"=>\$Help,
	   "shell:s"=>\$Shell,
	   "l:s"=>\$L,
          )
;
die `pod2text $0` if( $Help);
$L||=" vf=1G ";
my $comder="qsub   -cwd -S /bin/sh   -q ngb.q   -P ngb_un   -l " ;
$comder.=$L;
$comder.="   shell";


if($Shell){
          $comder.=" $Shell ";
	  }
else{	  
     if($ARGV[0]=~/perl/){
                         $comder.=" no  ";
		         }
     elsif($ARGV[0]=~/^sh/){
                             $comder.=" no  ";
                            }
     elsif($ARGV[0]=~/\.sh$/){
                              $comder.=" yes  "
	        	     }
     else{ die "commander wrong.";}	
     }


foreach(@ARGV){
              $comder.=" $_";
	      }
print "happy everyday:\n$comder\n";

`$comder`;


