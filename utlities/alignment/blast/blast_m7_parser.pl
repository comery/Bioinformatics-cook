#!/usr/bin/perl -w
use strict;

die "perl $0 <m7_file> <outfile>" unless @ARGV == 2;
my ($m7_file, $outfie) = @ARGV;

my ($quary_id, $quary_length, $query_start, $query_end, $subject_id, $subject_length, $subject_start, $subject_end, $identity, $positive, $align_length, $score, $evalue, $subject_annotation) = ""x14;
my $gap = 0;
open M7,"$m7_file";
open OUT,">$outfie";
my ($strand, $from, $to);
while(<M7>){
	chomp;
	if($_ =~ /<Iteration_query-def>(.*)<\/Iteration_query-def>/){
		$quary_id = $1;
	}
	elsif($_ =~ /<Iteration_query-len>(.*)<\/Iteration_query-len>/){
		$quary_length = $1;
	}
	elsif($_ =~ /<Hit_def>(.*)<\/Hit_def>/){
		my $tmp_str = &trim($1);
		if($tmp_str =~ /^(\S+)\s+(.*)$/){
			$subject_id = $1;
			$subject_annotation = $2;
			$subject_annotation =~ s/&gt;/>/g;
		}
	}
	elsif($_ =~ /<Hit_len>(.*)<\/Hit_len>/){
		$subject_length = $1;
	}
	elsif($_ =~ /<Hsp_bit-score>(.*)<\/Hsp_bit-score>/){
		$score = $1;
		if($score =~ /\./){
			$score = sprintf("%.1f", $score);;
		}
	}
	elsif($_ =~ /<Hsp_evalue>(.*)<\/Hsp_evalue>/){
		$evalue = $1;
		if($evalue =~ /(.*)e(.*)/){
			$evalue = int($1)."e".$2;;
		}
		else{
			$evalue = sprintf("%.1f", $evalue);
		}
	}
	elsif($_ =~ /<Hsp_query-from>(.*)<\/Hsp_query-from>/){
		$from = $1;
	}
	elsif($_ =~ /<Hsp_query-to>(.*)<\/Hsp_query-to>/){
		$to = $1;
	}
	elsif($_ =~ /<Hsp_hit-from>(.*)<\/Hsp_hit-from>/){
		$subject_start = $1;
	}
	elsif($_ =~ /<Hsp_hit-to>(.*)<\/Hsp_hit-to>/){
		$subject_end = $1;
	}
	elsif($_ =~ /<Hsp_query-frame>(.*)<\/Hsp_query-frame>/){
		$strand = ($1 =~ /-/) ? "-" : "+";
		if($strand eq "+"){
			$query_start = $from;
			$query_end = $to;
		}
		elsif($strand eq "-"){
			$query_start = $to;
			$query_end = $from;
		}
	}
	elsif($_ =~ /<Hsp_identity>(.*)<\/Hsp_identity>/){
		$identity = $1;
	}
	elsif($_ =~ /<Hsp_positive>(.*)<\/Hsp_positive>/){
		$positive = $1;
	}
	elsif($_ =~ /<Hsp_gaps>(.*)<\/Hsp_gaps>/){
		$gap = $1;
	}
	elsif($_ =~ /<Hsp_align-len>(.*)<\/Hsp_align-len>/){
		$align_length = $1;
		$identity = sprintf("%.2f", ($identity / $align_length));
		$positive = sprintf("%.2f", ($positive / $align_length));
#		print "$gap\t$align_length\n";
		$gap = sprintf("%.2f", ($gap / $align_length));
	}
	elsif($_ =~ /<\/Hit>/){
		print OUT "$quary_id\t$quary_length\t$query_start\t$query_end\t$subject_id\t$subject_length\t$subject_start\t$subject_end\t$identity\t$positive\t$gap\t$align_length\t$score\t$evalue\t--\t$subject_annotation\n";
		($query_start, $query_end, $subject_id, $subject_length, $subject_start, $subject_end, $identity, $positive, $align_length, $score, $evalue, $subject_annotation) = ""x12;
		$gap = 0;
	}
}
close M7;
close OUT;

sub trim{
	my ($in_str) = @_;
	$in_str =~ s/^\s+//;
	$in_str =~ s/\s+$//;
	return $in_str;
}
