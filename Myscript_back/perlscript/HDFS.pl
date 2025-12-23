#!/usr/bin/perl

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
  --in    <string>           the input file name           
  --out   <string>           the output file name
  --help                     output help information to screen  

=head1 Usage Exmples
  perl  ./syn_and_nsyn.pl  -start  1 -end 2280 -in Influenza-A-seg1 
  
=cut

#################################################################################
#use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max);

my ($start,$end,$Help,$cds_len,$inputfile,$outputfile);
GetOptions(
	"start:n"=>\$start,
	"end:n"=>\$end,
	"help"=>\$Help,
	"in:s"=>\$inputfile,
	"out:s"=>\$outputfile,
	) || die "Please use --help option to get help\n";
die `pod2text $0` if ($Help );

