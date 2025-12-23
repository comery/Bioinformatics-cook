#!/usr/bin/perl -w 
##############################################################
# This script is to change the format from genebank to fasta,
#precisely for CDS. It can split the joined CDS regions.
#Auther:
#		Chentao YANG yangchento@genomics.cn 2015/04/10
#############################################################
use strict;
use Bio::SeqIO;
use Bio::SeqFeatureI;
my $seq_obj;
open ERR,">err.log";
my $in = Bio::SeqIO-> new(-file => shift, "-format" => 'genbank');
while (my $seq_obj=$in->next_seq()) {
	my $primary_id = $seq_obj->primary_id;
	my $key = $seq_obj->accession_number;
	my $version = $seq_obj->version;
#	print "$primary_id";
#	my $gi = $1 if ($version =~ /GI:\s(\d+)/);

	for my $feature ($seq_obj->top_SeqFeatures){
		my ($val,$location,$subseq,$start,$end);
		if ( $feature->location->isa('Bio::Location::SplitLocationI')){
			if ($feature->primary_tag eq 'CDS' )  {
				if ($feature->has_tag('gene')) {
		   			for $val ($feature->get_tag_values('gene')){
		        #		print $feature->spliced_seq->seq,"\n";
						my @loca_count = ($feature->location->sub_Location);
						my $i = 1;
						my %hash = map{$_ => $i++}@loca_count;
				#		print "$loca_count\n";
						for $location ( $feature->location->sub_Location) {
							$start = $location->start;$end = $location->end;
		    				$subseq = $seq_obj->subseq($start,$end);
							my $j = $hash{$location};
							my $pro_name = $val."-$j";
		        			print ">gi|$primary_id|\t$pro_name\t$start\t$end\n";
							print $subseq."\n";
						}
					}
				}else{
					print ERR "Be careful! This feature don't have gene tag, you are supposed to add by hand!\n";
				}
			}
		}else{
			if ($feature->primary_tag eq 'CDS' ) {
				if ($feature->has_tag('gene')) {
		   			for $val ($feature->get_tag_values('gene')){
		        #		print $feature->spliced_seq->seq,"\n";
						for $location ( $feature->location) {
							$start = $location->start;$end = $location->end;
		    				$subseq = $seq_obj->subseq($start,$end);
		        			print ">gi|$primary_id|\t$val\t$start\t$end\n";
							print $subseq."\n";
						}
					}
				}else{
					print ERR "Be careful! This feature don't have gene tag, you are supposed to add by hand!\n";
				}
			}
		}
	}

}
close ERR;
