#!/usr/bin/perl -w

#########################################
#
# fasta2pfam: change fasta sequence format
#
# Usage:  aligned_TO_unaligned.pl <input.fasta> <output.pfam>
#
# Author: Chentao Yang yangchentao@genomics.cn  
# Date: March 24, 2015
##
########################################
use strict;
use Bio::Seq;
use Bio::SeqIO;
use IO::String;

#######################################
#
# Set up usage statement
#
#######################################
my $usage = " Usage:  aligned_TO_unaligned.pl <input.fasta> <output.pfam>\n";

my $scripthelp = "aligned_TO_unaligned.pl - change the aligned fasta sequence to be unaligned\n";

#######################################
#
# Initialize variables
#
#######################################
my $inFilename;
my $outFilename;
my $argNum = 2;

#######################################
#
# Test for commandline arguments
#
#######################################

if (! $ARGV[0] ) {
	print $usage;
	print $scripthelp;
	exit -1;
} elsif (scalar @ARGV != $argNum) {
	print $usage;
	print $scripthelp;
	exit -1;
} 

$inFilename = $ARGV[0];
if (! -f $inFilename) {
	print "Unable to locate input fasta file: $inFilename.\n";
	exit;
}
$outFilename = $ARGV[1];

