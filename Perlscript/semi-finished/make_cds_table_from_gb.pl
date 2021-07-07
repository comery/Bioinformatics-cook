use Bio::SeqIO;
use Bio::SeqFeatureI;
use Bio::DB::Taxonomy;
my $inputfile = shift or die "Usage: perl $0 <*.mito.gb>";
#open OUT ,">same_strand.log";
open KAR,">karyotype.txt";
my $in = Bio::SeqIO-> new(-file => "$inputfile", "-format" => 'genbank');
while (my $seq_obj=$in->next_seq()) {
	my $taxonomy = join(" ", $seq_obj->species->classification); 
	my $organism = $seq_obj->species->node_name; 
	$organism =~ s/ /\_/g;
	my $source = $seq_obj->seq;
	my $sou_len = $seq_obj->length;
	my $len_deN = $sou_len - $Ns;
	my $locus = $seq_obj -> display_id;

	print KAR "chr1 - mt1	mt1 0 $sou_len blue";

#	my @tmp;
	open CDS,">mitogenome_cds.txt";
	open RRNA,">mitogenome_rrna.txt";
	open TRNA,">mitogenome_trna.txt";
	open TEXT,">mitogenome_text.txt";

	for my $feature ($seq_obj->top_SeqFeatures){
		my ($db_xref,$val,$location,$subseq,$pro,$translation);

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
				$direction = '?';
			}
			push @tmp,$direction;
			print CDS "mt1\t$start\t$end\n";

			if ($feature->has_tag('gene')) {
				for $val ($feature->get_tag_values('gene')){
				
					print TEXT "mt1\t$start\t$end\t$val\n";
				}
			}
		}elsif ($feature->primary_tag eq 'rRNA' )  {
		
			my $start = $feature -> start;
			my $end = $feature -> end;
			print RRNA "mt1\t$start\t$end\n";

			if ($feature->has_tag('gene')) {
				for $val ($feature->get_tag_values('gene')){
				
					print TEXT "mt1\t$start\t$end\t$val\n";
				}
			}else{
			
				print TEXT "mt1\t$start\t$end\tNA\n";
			}

		}elsif ($feature->primary_tag eq 'tRNA' )  {
		 
		 	my $start = $feature -> start;
			my $end = $feature -> end;
			print TRNA "mt1\t$start\t$end\n";
			if ($feature->has_tag('gene')) {
				for $val ($feature->get_tag_values('gene')){
				
					print TEXT "mt1\t$start\t$end\t$val\n";
				}
			}else{
				
					print TEXT "mt1\t$start\t$end\tNA\n";
			}
		}
	}
	
}

close KAR;
close CDS;
close RRNA;
close TRNA;
close TEXT;

# creat a circos.conf file
open CON,">circos.conf";
print  CON <<_CONFIG_;

<<include etc/colors_fonts_patterns.conf>>

<<include ideogram.conf>>
<<include ticks.conf>>

show_ticks* = no

karyotype   = karyotype.txt

<image>
<<include etc/image.conf>>
</image>



chromosomes_units = 1000000
chromosomes       = mt1
chromosomes_display_default = no

<plots>

type       = text
color      = black
label_font = default
label_size = 24p

#GC content
<plots>
<<include plots_histogram.conf>>  
</plots>

#description info
<plot>
file = mitogenome_text.txt

r1   = 1r+200p
r0   = 1r

show_links     = yes
link_dims      = 0p,0p,70p,0p,10p
link_thickness = 2p
link_color     = red

</plot>


</plots>

<highlights>

# CDS
<highlight>
file         = mitogenome_cds.txt
r0           = 1.09r
r1           = 1.11r
fill_color   = orange
z            = 5
</highlight>

# rRNA
<highlight>
file         = mitogenome_rrna.txt
r0           = 0.90r
r1           = 0.94r
fill_color   = green
z            = 4
stroke_thickness = 1p
</highlight>

# tRNA
<highlight>
file         = mitogenome_trna.txt
r0           = 0.90r
r1           = 0.94r
fill_color   = purple
z            = 2
stroke_thickness = 1p
</highlight>

s = 1
</highlight>


</highlights>


<<include etc/housekeeping.conf>>

_CONFIG_

close CON;
