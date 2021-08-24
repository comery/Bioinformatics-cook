#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);
use Bio::SeqIO;
use Bio::SeqFeatureI;
use Getopt::Long;
use Statistics::Descriptive;

=head1 Description
	This script is an utility for chloroplast/mitochondria 
	annotation pipoline to draw circos pictures 

=head1 Author :
	YangChentao yangchentao@genomics.cn 2017/05

=head1 options
	*--gb 	<str>	*.mitogenome.gb
	--outdir <gc>	./gc
	--help		display this help information

=head1 Usage
	perl draw_circos_for_mitogenome.pl  -gb <*.mitogenome.gb>  -outdir gc

=head1 attention
	

=cut


my ($help,$gbfile,$warn,$outdir);

GetOptions(
			'help' => \$help,
			'gb=s' => \$gbfile,
			'outdir=s' => \$outdir
	);

die `pod2text $0` if ($help || !$gbfile );

$outdir ||= 'Class_GC';
$outdir = abs_path($outdir);


`[ -d $outdir ] || mkdir $outdir`;
print "outdir : $outdir...\n";
### read genebank file 
# whether link chromsome as a whole according to topology[linear|circular] of mitogenome.

# get cds rRNA tRNA's location and strand information
my $in = Bio::SeqIO-> new(-file => "$gbfile", "-format" => 'genbank');
while (my $seq_obj=$in->next_seq()) {
	
	my $source = $seq_obj->seq;
	my $sou_len = $seq_obj->length;
	my $locus = $seq_obj -> display_id;

	# class name is prefix of $locus, making hash of groups
	my $class = (split /_/,$locus)[0];
	`mkdir -p $outdir/$class` unless ( -d "$outdir/$class");

#	record the class groups information
	
	#calculate GC content
#	if ($conf{'gc'} eq 'yes') {
	open GC,">$outdir/$class/$locus.fa.gc.table";

	for my $feature ($seq_obj->top_SeqFeatures){
		my ($db_xref,$val,$location,$gc_num,$gc,$seq,$start,$end,$G,$C,$gg,$cc);

		if ($feature->primary_tag eq 'CDS' )  {
			$seq = $feature->spliced_seq->seq;
			$start = $feature -> start;
			$end = $feature -> end;
			$gc_num = $seq =~ tr/GCgc/GCgc/;
			$G = $seq =~ tr/Gg/Gg/;
			$C = $seq =~ tr/Cc/Cc/;
			$gc = $gc_num/(abs($end - $start) + 1);
			$gc = sprintf("%.4f",$gc);
			$gg = $G/$gc_num;$gg = sprintf("%.4f",$gg);
			$cc = $C/$gc_num;$cc = sprintf("%.4f",$cc);
			if ($feature->has_tag('gene')) {
				for $val ($feature->get_tag_values('gene')){
				
					print GC "$val\t$gc\t$gg\t$cc\n";
				}
			}elsif ($feature->has_tag('product')){
				for $val ($feature->get_tag_values('product')){
				
					print GC "$val\t$gc\t$gg\t$cc\n";
				}
			}else {
				print GC "CDS_NA\t$gc\t$gg\t$cc\n";
			}
		}elsif ($feature->primary_tag eq 'rRNA' )  {
		
			$seq = $feature->spliced_seq->seq;
			$start = $feature -> start;
			$end = $feature -> end;
			$gc_num = $seq =~ tr/GCgc/GCgc/;
			$G = $seq =~ tr/Gg/Gg/;
			$C = $seq =~ tr/Cc/Cc/;
			$gc = $gc_num/(abs($end - $start) + 1);
			$gc = sprintf("%.4f",$gc);
			$gg = $G/$gc_num;$gg = sprintf("%.4f",$gg);
			$cc = $C/$gc_num;$cc = sprintf("%.4f",$cc);
		#	print RRNA "mt1\t$start\t$end\n";

			if ($feature->has_tag('gene')) {
				for $val ($feature->get_tag_values('gene')){
					print GC "$val\t$gc\t$gg\t$cc\n";
				}
			}elsif($feature->has_tag('product')){
				for $val ($feature->get_tag_values('product')){
					$val =~ s/\s/_/g;
					print GC "$val\t$gc\t$gg\t$cc\n";
				}
				
			}else {
				print GC "rRNA_NA\t$gc\t$gg\t$cc\n";
			}

	#	}elsif ($feature->primary_tag eq 'tRNA' )  {
	#	 	$seq = $feature->spliced_seq->seq;
	#		$start = $feature -> start;
	#		$end = $feature -> end;
	#		$gc_num = $seq =~ tr/GCgc/GCgc/;
	#		$gc = $gc_num/(abs($end - $start) + 1);
	#		$gc = sprintf("%.4f",$gc);
			
	#		if ($feature->has_tag('gene')) {
	#			for $val ($feature->get_tag_values('gene')){
				
	#				print GC "$val\t$gc\n";
	#			}
	#		}elsif($feature->has_tag('product')){
	#			for $val ($feature->get_tag_values('product')){
	#				print GC "mt1\t$start\t$end\t$val\n";
	#			}
	#		}else{
				
	#				print GC "tRNA_NA\t$gc\n";
	#		}
		}
	}
	

close GC;

}


print "Reading genebank file done. \n";
