#!/usr/bin/perl -w
use strict;
use Bio::DB::Taxonomy;
use Bio::DB::Taxonomy::flatfile;
# use NCBI Entrez over HTTP
my $usage = "perl $0 <species list>[Each line is a specimen like \"Homo sapiens\" ]";
die "$usage" unless (@ARGV > 0);
open IN,shift;
while (<IN>) {
	chomp;
	my $db = Bio::DB::Taxonomy->new(-source => 'flatfile',
									-nodesfile => '/ifs4/NGB_ENV/USER/yangchentao/software/perlscript/taxdmp/nodes.dmp',
	                                -namesfile => '/ifs4/NGB_ENV/USER/yangchentao/software/perlscript/taxdmp/names.dmp');
	my $taxonid = $db->get_taxonid('$_');
	print "$taxonid\n";
    #  # get a taxon
    my $taxon = $db->get_taxon(-taxonid => $taxonid);
    print "$_\t$taxon\n";
}

close IN;
