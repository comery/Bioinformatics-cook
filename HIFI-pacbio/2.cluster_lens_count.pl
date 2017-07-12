#!/usr/bin/perl -w
=head1 Name
	
	cluster_lens_count.pl -cluster ccs sequence by length and count and codon check

=head1 Usage
	
	perl cluster_lens_count.pl  [option]
	
	--ccs <str> all ccs fasta
	-pattern <str> check.ccs.fa.log
	--min <number> the min cutoff length for clustering
	--max <number> the max cutoff length for clustering
	--codon <number> the genetic codon table id,default 5(Invertebrate Mitochondrial)
	--frame <number> '0|1|2' ORF start shift,default 1

=head1 Exmple

	perl cluster_lens_count.pl -ccs ccs.fa.out -pattern check.ccs.fa.log -min 652 -max 664

=cut

use strict;
use Getopt::Long;
my ($ccs,$pattern,$passes,$minL,$maxL,$codon,$frame,$help);
GetOptions (
		"ccs=s" => \$ccs,
		"pattern=s" => \$pattern,
		"passes=s" => \$passes,
		"min:i" =>\$minL,
		"max:i" =>\$maxL,
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
	$/="//";
	my $str = <PA>;
	chomp $str;
	$str = (split /\n/,$str)[1];
	if ($a[1] == 1 && $str =~ /^seq_location:\d+-(\d+)/){$check{$a[0]}{'s'} = $1;}
	if ($a[1] == 2 && $str =~ /^seq_location:\d+-(\d+)/){$check{$a[0]}{'s'} = $1;}
	if ($a[1] == 3 && $str =~ /^seq_location:(\d+)-\d+/){$check{$a[0]}{'e'} = $1;}
	if ($a[1] == 4 && $str =~ /^seq_location:(\d+)-\d+/){$check{$a[0]}{'e'} = $1;}
	$/="\n";
}
close PA;
###---------------------------------------------------------------------------------

open LOG, ">cluster.all.fa";
open OUT , ">cluster.top1.fas";
open OUT1, ">cluster.top3.fas";

my (%cluster,%count,%seqs,$err_len,$allseq,$err_codon);
### trim primer region and keep target sequence(most are 658bp)
open CCS,"$ccs" or die "Can not open $ccs!";
$/=">";<CCS>;$/="\n";
while (my $cid=<CCS>) {
	chomp $cid;
#	print "$cid\n";
	my @a = split /\s+/,$cid;
	my $ccs_id = $a[0];
	my $sam = $a[1];
	my $ori = $a[2];
	$/=">";
	my $seq = <CCS>;
	chomp $seq;
	$seq =~ s/\n//g;
	$seq =~ s/\s//g;
#	print "$cid_id\n$seq\n";
	$seq = &trim_primer($seq,$check{$ccs_id}{'s'},$check{$ccs_id}{'e'});
	$seqs{$ccs_id} = $seq;
	my $l = length $seq;

	$allseq++;
	if ($l>=$minL && $l<=$maxL){
		
		if (&TranslateDNASeq($seq) >= 0) {
		#	print &TranslateDNASeq($seq)."\n";
		#	print "$cid has stop condon!\n";
			$err_codon ++;
		}else {
			$count{$sam} ++;
			if ($ori eq "for") {
				$cluster{$sam}{$seq} ++;
			}else {
				my $rev = &comrev($seq);
				$cluster{$sam}{$rev} ++;
			}
		}
		
	}else {
		$err_len++;
	}
	$/="\n";

}

print "$allseq ccs\n$err_len has wrong length\n$err_codon has wrong condon\n";

### reading passes number info-----
open PS, "$passes";
my %pass;
while (<PS>) {
	chomp;
	my @ps = split;
#	print "$seqs{$ps[0]}\n";
	$pass{$seqs{$ps[0]}} = $ps[1] if (exists $seqs{$ps[0]});
}
close PS;
###--------------------------------

my @sample = sort{$a <=> $b} keys %cluster;
foreach my $sam(@sample) {
	my $sort = $cluster{$sam};
	my @sorted = sort {$$sort{$b} <=> $$sort{$a} } keys %$sort;
	my $sam_num = @sorted;
	### sort by passes number if each kind of sequence's count is 1
	if ($cluster{$sam}{$sorted[0]} == 1) {
		@sorted = sort {$pass{$b} <=> $pass{$a}} @sorted;

	}
		#print best three sequences
		print OUT ">$sam\_1\n$sorted[0]\n";
		print OUT1 ">$sam\_1\n$sorted[0]\n" ;
		print OUT1 ">$sam\_2\n$sorted[1]\n" if ($sorted[1]);
		print OUT1 ">$sam\_3\n$sorted[2]\n" if ($sorted[2]);
	
		my $i = 0;
		foreach my $iso (@sorted) {
			my $len = length $iso;
			my $j = $i + 1;
			print LOG ">$sam\_$j\t$cluster{$sam}{$iso}/$count{$sam}\t$len\t$pass{$iso}\n$iso\n";
			$i++;
		}		
}


close CCS;
close OUT;
close OUT1;
close LOG;



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
