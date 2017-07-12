#!/usr/bin/perl -w
=head1 Name
	
	cluster_lens_count.pl -cluster ccs sequence by length and count and codon check

=head1 Usage
	
	perl cluster_lens_count.pl  [option]
	
	--ccs <str> all ccs fasta
	-pattern <str> check.ccs.fa.log
	--min <number> the min cutoff length for clustering
	--max <number> the max cutoff length for clustering
	--c whether to chech codon
	--codon <number> the genetic codon table id,default 5(Invertebrate Mitochondrial)
	--frame <number> '0|1|2' ORF start shift,default 1

=head1 Exmple

	perl cluster_lens_count.pl -ccs ccs.fa.out -pattern check.ccs.fa.log -min 652 -max 664 -c

=cut

use strict;
use Getopt::Long;
my ($ccs,$pattern,$passes,$minL,$maxL,$check,$codon,$frame,$help);
GetOptions (
		"ccs=s" => \$ccs,
		"pattern=s" => \$pattern,
		"passes=s" => \$passes,
		"min:i" =>\$minL,
		"max:i" =>\$maxL,
		"c!" => \$check,
		"codon:i" =>\$codon,
		"frame:i" =>\$frame,
		"help" =>\$help
);
die `pod2text $0` if ( !$ccs || ! $pattern || !$passes ||$help);
die "translation frame is wrong, it must be (0|1|2)" if ( $frame && $frame != 0 && $frame != 1 && $frame != 2);
$minL ||= 652; # 658 - 6bp
$maxL ||= 664; # 658 + 6bp
$codon ||= 5; # Invertebrate Mitochondrial
$frame ||= 1; # COI 658bp barcode from second base to translate, now I didn't consider that condition indel occured before the start codon

## collecting primer location info------------------------------------------------
open PA, "$pattern";
my %check;
$/="//";<PA>;$/="\n";
while (<PA>) {
	next if (/^#/);
	chomp;
	my @a = split;
	my $newid = "ccs_".(split /\//,$a[0])[-2];
	$/="//";
	my $str = <PA>;
	chomp $str;
	$str = (split /\n/,$str)[1];
	if ($a[1] == 1 && $str =~ /^seq_location:\d+-(\d+)/){$check{$newid}{'s'} = $1;}
	if ($a[1] == 2 && $str =~ /^seq_location:\d+-(\d+)/){$check{$newid}{'s'} = $1;}
	if ($a[1] == 3 && $str =~ /^seq_location:(\d+)-\d+/){$check{$newid}{'e'} = $1;}
	if ($a[1] == 4 && $str =~ /^seq_location:(\d+)-\d+/){$check{$newid}{'e'} = $1;}
	$/="\n";
}
close PA;
###---------------------------------------------------------------------------------

open ALL, ">cluster.all.fa";
open OUT , ">cluster.top1.fas";
open LOG, ">cluster.id.txt";

print LOG "sample_top\tcount/all_ccs_number\tlength\ttotal_passes_number\tccs_id_cluster\n";

my (%cluster,%cluster_id,%count,%seqs,$err_len,$allseq,$err_codon);
### trim primer region and keep target sequence(most are 658bp)
open CCS,"$ccs" or die "Can not open $ccs!";
$/=">";<CCS>;$/="\n";
while (my $cid=<CCS>) {
	chomp $cid;
	my @a = split /\s+/,$cid;
	my $ccs_id = "ccs_".(split /\//,$a[0])[-2];
	my $sam = $a[1];
	my $ori = $a[2];
	$/=">";
	my $seq = <CCS>;
	chomp $seq;
	$seq =~ s/\n//g;
	$seq =~ s/\s//g;
	$seq = &trim_primer($seq,$check{$ccs_id}{'s'},$check{$ccs_id}{'e'});
	$seqs{$ccs_id} = $seq;
	my $l = length $seq;

	$allseq++;

	my $judge = (defined $check ? (&TranslateDNASeq($seq) >= 0) : 0 );

	if ($l>=$minL && $l<=$maxL){
		
		if ($judge) {
		#	print &TranslateDNASeq($seq)."\n";
		#	print "$cid has stop condon!\n";
			$err_codon ++ if (defined $check);
		}else {
			$count{$sam} ++;
			if ($ori eq "for") {
				$cluster{$sam}{$seq} ++;
				push @{$cluster_id{$seq}},$ccs_id;
			}else {
				my $rev = &comrev($seq);
				$cluster{$sam}{$rev} ++;
				push @{$cluster_id{$rev}},$ccs_id;
			}
		}
		
	}else {
		$err_len++;
	}
	$/="\n";

}

if (defined $check) {
	print "chech codon model!\n";
	print "total ccs:$allseq\n$err_len has wrong length\n$err_codon has wrong condon\n" ;
}else {
	print "do not check codon model!\n";
	print "total ccs:$allseq\n$err_len has wrong length\n";
}

### reading passes number info-----
open PS, "$passes";
my %pass;
while (<PS>) {
	chomp;
	my @ps = split;
	my $id = "ccs_".(split /\//,$ps[0])[-2];
	$pass{$id} = $ps[1];
}
close PS;
###--------------------------------

### statistic total passes number of each kind of CCS
my %passes_dis;
foreach (keys %cluster_id) {
	my $total = &total_pass(@{$cluster_id{$_}}) ;
	$passes_dis{$_} = $total;
}
###-------------------------------------------------

### print out
my @sample = sort{$a <=> $b} keys %cluster;
foreach my $sam(@sample) {
	my $sort = $cluster{$sam};
	### sort by total passes number
	my @sorted = sort {$passes_dis{$b} <=> $passes_dis{$a}} keys %$sort;
		#print best three sequences
		my $top_len = length $sorted[0];
		print OUT ">$sam\_1\tcount=$cluster{$sam}{$sorted[0]}/$count{$sam}\tlen=$top_len\ttotal_passes=$passes_dis{$sorted[0]}\n$sorted[0]\n";
	
		my $i = 0;
		foreach my $iso (@sorted) {
			my $len = length $iso;
			my $j = $i + 1;
			print ALL ">$sam\_$j\tcount=$cluster{$sam}{$iso}/$count{$sam}\tlen=$len\ttotal_passes=$passes_dis{$iso}\n$iso\n";
			print LOG "$sam\_$j\t$cluster{$sam}{$iso}/$count{$sam}\t$len\t$passes_dis{$iso}\t[@{$cluster_id{$iso}}]\n";
			$i++;
		}		
}


close CCS;
close OUT;
close LOG;
close ALL;



sub comrev {
	my $a = shift;
	chomp $a;
	$a =~ tr/NATCG/NTAGC/;
	$a= reverse $a;
	return $a;
}

sub trim_primer {
	my @a = @_;
	my $seq = $a[0];
	my $s = $a[1];
	my $e = $a[2];
	my $len = length $seq;
	my $end_trim = $len - $e +1;
	substr($seq,0,$s)  = "";
	substr($seq,-$end_trim) = "";
	return $seq;
}


sub TranslateDNASeq()
   {
      use Bio::Seq;
      (my $dna)=@_;
      my $seqobj=Bio::Seq->new(-seq =>$dna, -alphabet =>'dna');
      my $prot_obj = $seqobj->translate(-codontable_id => $codon,-terminator => 'U',-unknown => '_',-frame => $frame);
      my $pep = $prot_obj -> seq();
      my $stop_first = index($pep,"U");
      return $stop_first;
   }

sub total_pass {

	my $t;
	foreach (@_) {
		$t += $pass{$_};
	}
	return $t;
}
