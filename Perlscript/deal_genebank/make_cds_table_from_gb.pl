use Bio::SeqIO;
use Bio::SeqFeatureI;
use Bio::DB::Taxonomy;
my $inputfile = shift or die "Usage: perl $0 <*.mito.gb>";
#open OUT ,">same_strand.log";
my $in = Bio::SeqIO-> new(-file => "$inputfile", "-format" => 'genbank');
while (my $seq_obj=$in->next_seq()) {
	my $keywd = $seq_obj->desc();
	print "$keywd\n";
	my $taxonomy = join(" ", $seq_obj->species->classification); 
	my $organism = $seq_obj->species->node_name; 
	$organism =~ s/ /\_/g;
	my $source = $seq_obj->seq;
	my $sou_len = $seq_obj->length;
	my $Ns;
	$source =~ /N/ ? $Ns = $source =~ s/N/N/g : $Ns = 0;
	my $len_deN = $sou_len - $Ns;
	my $locus = $seq_obj -> display_id;
#	my @tmp;
	for my $feature ($seq_obj->top_SeqFeatures){
		my ($db_xref,$val,$location,$subseq,$start,$end,$pro,$translation);

		if ($feature->primary_tag eq 'CDS' )  {
			my $seq = $feature->spliced_seq->seq;
			my $start = $feature -> start;
			my $end = $feature -> end;
			my $strand = $feature -> strand;
			my $direction;
			if ($strand == 1) {
				$direction = '+';
			}elsif ($strand == -1){
				$direction = '-';
			}else {
				$direction = 'no relevant';
			}
			push @tmp,$direction;
		#           $feature->has_tag('gene') ? $val = ($feature->get_tag_values('gene'))[0] : $val = "NA";
		#           $feature->has_tag('product') ? $pro = ($feature->get_tag_values('product'))[0] : $pro = "NA";
			if ($feature->has_tag('gene')) {
				for $val ($feature->get_tag_values('gene')){
					print "$organism\t$val\t$locus\t$sou_len\t$len_deN\t$start\t$end\t$direction\n";
				}
			}
		}
	}
	#print OUT "$organism\t@tmp\n";
	#my %count;
#	my @uniq_tmp = grep {++$count{$_} < 2} @tmp;
	#print OUT "$organism\tall genes are in $uniq_tmp[0] strand!\n" if (@uniq_tmp = 1);
	#@tmp = ();
}
