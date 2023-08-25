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

#######################################
#
# Open the files
#
#######################################

my $in = Bio::SeqIO->new('-file' => "<$inFilename", '-format' => "fasta") || die ("Unable to read input fasta file: $inFilename.  Exiting.\n");

open (OUT, ">$outFilename") || die "Unable to write to $outFilename\n";
#print OUT "# STOCKHOLM 1.0\n";

while (my $seqobj = $in->next_seq)
{
	my $str = $seqobj->seq;
	$str =~ s/[\_-]//g;
	my $id = $seqobj->id;
	print OUT ">$id\n$str\n";
}

#print OUT "//\n";
