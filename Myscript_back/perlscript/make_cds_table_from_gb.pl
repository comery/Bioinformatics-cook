
use Bio::SeqIO;
use Bio::SeqFeatureI;
my $inputfile = shift or die "$!";
my $in = Bio::SeqIO-> new(-file => "$inputfile", "-format" => 'genbank');
while (my $seq_obj=$in->next_seq()) {
	my $source = $seq_obj->seq;
	my $sou_len = $seq_obj->length;
	my $Ns;
	$source =~ /N/ ? $Ns = $source =~ s/N/N/g : $Ns = 0; 
	my $len_deN = $sou_len - $Ns;
	my $locus = $seq_obj -> display_id;
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
#           $feature->has_tag('gene') ? $val = ($feature->get_tag_values('gene'))[0] : $val = "NA";
#           $feature->has_tag('product') ? $pro = ($feature->get_tag_values('product'))[0] : $pro = "NA";

        if ($feature->has_tag('gene')) {
            for $val ($feature->get_tag_values('gene')){
                print "$val\t$locus\t$sou_len\t$len_deN\t$start\t$end\t$direction\n";
            }
        }
    }
}
}

           
