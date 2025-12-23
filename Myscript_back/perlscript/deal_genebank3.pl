#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::SeqFeatureI;
use Bio::Species;
my $seq_obj;
my $in = Bio::SeqIO-> new(-file => shift, "-format" => 'genbank');
while (my $seq_obj=$in->next_seq()) {
    my $id = $seq_obj->accession_number;
    my $seq = $seq_obj->seq;
    my $seq_length = $seq_obj->length;
    my $taxname = $seq_obj->species->node_name;
	$taxname =~ s/ /\_/g;
	`mkdir -p species/$taxname`;
	open TRNA,">>species/$taxname/tRNA.fa";
	open RRNA,">>species/$taxname/rRNA.fa";
	open CDS,">>species/$taxname/CDS.fa";
    
	for my $feat ($seq_obj->get_SeqFeatures) {
        if ($feat->primary_tag eq "rRNA") {
            my $r_subseq = $feat->spliced_seq->seq;
            my $r_length = length($r_subseq);
			my $r_val;
        	if ($feat->has_tag('gene')) {
				$r_val = ($feat->get_tag_values('gene'))[0];
			}else{
				$r_val = "NA";
			}
			print RRNA ">gi_".$id."_".$r_val."_".$taxname."_".$r_length."_rRNA\n".$r_subseq."\n";
		}elsif($feat->primary_tag eq "tRNA"){
            my $t_subseq = $feat->spliced_seq->seq;
            my $t_length = length($t_subseq);
			my $t_val;
            if ($feat->has_tag('gene')) {
				$t_val = ($feat->get_tag_values('gene'))[0];
			}else{
				$t_val = "NA";
			}
			print TRNA ">gi_".$id."_".$t_val."_".$taxname."_".$t_length."_tRNA\n".$t_subseq."\n";
        }elsif($feat->primary_tag eq "CDS"){
            my $cd_subseq = $feat->spliced_seq->seq;
            my $cd_length = length($cd_subseq);
			my ($val,$produce);
            if ($feat->has_tag('gene')) {
				$val = ($feat->get_tag_values('gene'))[0];
			}else{
				$val = "NA";
			}
            if ($feat->has_tag('product')) {
				$produce = ($feat->get_tag_values('product'))[0];
			}elsif($feat->has_tag('note')){
				$produce = ($feat->get_tag_values('note'))[0];
			}else{
				$produce = "NA";
			}
            print CDS ">gi_".$id."_".$val."_".$taxname."_".$cd_length."_"."$produce\n".$cd_subseq."\n";

        }
    }

close RRNA;
close TRNA;
close CDS;
}
#sub joinseq {
#	my $feature = shift;
#	my $join;
#	if ( $feature->location->isa('Bio::Location::SplitLocationI')){
#  		for my $location ( $feature->location->sub_Location ) {
#			$join .= $feature->subseq($location->start,$location->end);
#      }
#   }
#   return $join;
#}
