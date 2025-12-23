#!/usr/bin/perl -w
use strict;
use Bio::DB::Taxonomy;
use Bio::DB::Taxonomy::flatfile;
# use NCBI Entrez over HTTP
#my $usage = "perl $0 <species list>[Each line is a specimen like \"Homo sapiens\" ]";
#die "$usage" unless (@ARGV == 0);
open IN,shift;
my $obj = Bio::DB::Taxonomy::flatfile->new(-source => 'flatfile',
                                            -directory => 'tax_data',
                                            -force     => 1 ,
                                          -nodesfile => 'tax_data/nodes.dmp',
                                          -namesfile => 'tax_data/names.dmp');


while (<IN>) {
  chomp;
#  my $db = $obj;
  my $taxonid = $obj->get_taxonid('$_');
#  # get a taxon
  my $taxon = $obj->get_taxon(-taxonid => $taxonid);
  print "$_\t$taxon\n";
}

close IN;

