#!/usr/bin/perl -w 
use strict;
use Bio::SeqIO;
use Bio::SeqFeatureI;
use Bio::SeqFeature::Generic;
my $seq_obj;
die "Usage:\n\tperl $0 <source.fasta>  <most_annot.sorted.cds> <output>" unless (@ARGV == 3);
my ($source,$features,$output);
$source = $ARGV[0];
$features = $ARGV[1];
$output = $ARGV[2];
open FA, $source;
open FEA,$features;
## Make a hash for scaffolds ##
$/= ">";<FA>;$/="\n";
my %scaff;
while (my $id = <FA>){
	chomp $id;
	$id = (split /\s+/,$id)[0];
	$/=">";
	my $str = <FA>;
	chomp $str;
	$scaff{$id} = $str;
	$/="\n";
}

## Make a two-dimensional hash for features of scaffold ##
my %feats;
while (my $string = <FEA>) {
	my $location = (split /\s+/,$string)[2];
	my $feat_id = $1 ,my $feat_s = $2 if ($location =~ /locus=(\w+):(\d+):(\d+):/);
	$feats{$feat_id}{$feat_s} = $string;
}

foreach my $key (keys %scaff) {
	my $seq = $scaff{$key};
	my $len = length($seq);
	my $seq_obj = Bio::Seq->new(-seq => $seq,
								-locus => "$key\tbp\tDNA",
								-difinition => "Green plant chloroplast",
								-author => "Yang Chentao",
	                            -display_id => "$output" );
	my $feat_source = new Bio::SeqFeature::Generic(	-start	=>	1,
													-end	=>	$len,
													-primary_tag => 'source',
													-tag	=>	{-organism	=>	'Weed',
																 -organelle	=>	"plastid:chloroplast",
																 -mol_type	=>	'chloroplast genome'}
												   );
	$seq_obj->add_SeqFeature($feat_source);
	my $hash2 = $feats{$key};
	foreach my $key2 (keys %$hash2) {
		my $val = $$hash2{$key2};
		my @aa = split /\s+/,$val;
		my @info = split /_/,$aa[0];
		my @bb = split /:/,$aa[2];
		my $genename = $info[2];
		my $note = "$info[3] $info[4] $info[5]";
		my $start = $bb[1];
		my $end = $bb[2];
		my $strand;
		if ($bb[3] eq "+" ) {
			$strand = 1 ;
		}else{
			$strand = -1 ;
		}
		
		# create the feature with some data, evidence and a note
		my $feat = new Bio::SeqFeature::Generic(-start       => $start,
		                                        -end         => $end,
											    -strand      => $strand,
											    -primary_tag => 'CDS',
												-tag => {gene => $genename,note => $note } );
		$seq_obj->add_SeqFeature($feat);
	}
	my $io = Bio::SeqIO->new(-format => "genbank", -file => ">$output.gb" );       
	$io->write_seq($seq_obj);
}	
