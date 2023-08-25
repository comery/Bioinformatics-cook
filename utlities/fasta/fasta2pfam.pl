#!/usr/bin/perl -w

#########################################
#
# fasta2pfam: convert from fasta to pfam
#
# Usage:  fasta2pfam <input.fasta> <output.pfam>
#
# Author: Susan Huse, shuse@mbl.edu  
# Date: January 11, 2005
#
# Copyright (C) 2005 Marine Biological Laborotory, Woods Hole, MA
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# For a copy of the GNU General Public License, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# or visit http://www.gnu.org/copyleft/gpl.html
#
# Keywords: convert pfam
# 
# Assumptions: 
#
# Revisions:
#
# Programming Notes:
#
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
my $usage = " Usage:  fasta2pfam <input.fasta> <output.pfam>
                 ex:  fasta2pfam in.fasta out.pfam\n";
my $scripthelp = "fasta2pfam - converts a fasta file and pfam format\n";

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
print OUT "# STOCKHOLM 1.0\n";

while (my $seqobj = $in->next_seq)
{
	my $str = $seqobj->seq;
	my $id = $seqobj->id;
	print OUT "$id $str\n";
}

print OUT "//\n";
