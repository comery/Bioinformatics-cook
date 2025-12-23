#!/usr/local/bin/perl -w
use strict;
my $data = $ARGV[0];
#open IN,$data;
my ($in,$out);
use Bio::AlignIO;
my $inputfilename = "$data";

$in = Bio::AlignIO->new(-file => $inputfilename , '-format' => 'fasta');
$out = Bio::AlignIO->new(-file => ">$data.nex" , '-format' => 'nexus');
#format included "nexus fasta "
# note: we quote -format to keep older perl's from complaining.

while ( my $aln = $in->next_aln() ) {
  $out->write_aln($aln);
}

