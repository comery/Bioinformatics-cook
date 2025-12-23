use Bio::SearchIO;
my $blast_result = shift;
   # format can be 'fasta', 'blast', 'exonerate', ...
   my $searchio = Bio::SearchIO->new( -format => 'blastxml',
                                     -file   => $blast_result );
   while ( my $result = $searchio->next_result() ) {
     my $query_name = $result->query_name;
     my $query_len = $result->query_length;
     while( my $hit = $result->next_hit ) {
       my $hit_name = $hit->name;
       my $des = ($hit->description ? $hit->description : 'NA');
        # process the Bio::Search::Hit::HitI object
           while( my $hsp = $hit->next_hsp ) { 
             # process the Bio::Search::HSP::HSPI object
             my $e = $hsp->evalue;
             my $identity = $hsp->frac_identical;
			 my $aln_len = $hsp->hsp_length;
			 my $conserved = $hsp->num_conserved;
			 my $identical = $hsp->num_identical;
			 my $mismatch = $hsp->length('hit') - $identical; ## length of hit participating in alignment minus gaps,then minus identical base;
			 my $gap = $hsp->gaps;
			 my $score = $hsp->score;
             my $start_q = $hsp->start('query');
             my $end_q = $hsp->end('query');
             my $start_h = $hsp->start('hit');
             my $end_h = $hsp->end('hit');
			 use Bio::AlignIO;

			 # $aln will be a Bio::SimpleAlign object
			 my $aln = $hsp->get_aln; 
			 my $alnIO = Bio::AlignIO->new(-format => "msf", 
			                                -file => ">hsp.msf"); 
			$alnIO->write_aln($aln);
             print "$query_name\t$hit_name\t$identity\t$aln_len\t$mismatch\t$gap\t$start_q\t$end_q\t$start_h\t$end_h\t$e\t$score\t$des\n";

            
           }
       }
   }
