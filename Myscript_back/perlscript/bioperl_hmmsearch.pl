
use strict;
use Bio::SearchIO;
my $in = new Bio::SearchIO(-format => 'hmmer',-file => shift);
my ($hit,$hsp);
while ( my $res = $in->next_result ){
	# get a Bio::Search::Result::HMMERResult object
	print $res->query_name, " for HMM ", $res->hmm_name, "\n";
	while ( $hit = $res->next_hit ){
		print $hit->name, "\n";
		while ( $hsp = $hit->next_hsp ){
			print "length is ", $hsp->length, "\n";
		}
	}
}
