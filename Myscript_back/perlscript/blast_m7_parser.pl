#!/usr/bin/perl -w

=pod
description: convert blast's result from -m7 (xml) to -m8 (tabular, without comment lines)
author: Zhang Fangxian, zhangfx@genomics.cn
created date: 20090910
modified date: 20101127, 20091029, 20091025
=cut

use strict;

my ($input, $output) = @ARGV;

open OUT, "> $output" or die $!;
print OUT "Query_id\tSubject_id\tIdentity\tAlign_length\tMiss_match\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE_value\tScore\tSubject_annotation\n";
open IN, $input or die $!;
	my $gap = 0;
	my ($query, $hitId, $hitDef, $bits, $evalue, $qFrom, $qTo, $hFrom, $hTo, $frame, $identity, $length, $annot);
	while (<IN>) {
		if (/<(Iteration_query-def)>(\S+).*<\/\1/) {
			$query = $2;
		} elsif (/<(Hit_id)>(.*)<\/\1/) {
			$hitId = $2;
		} elsif (/<(Hit_def)>([^\s]*)\s+(.*)<\/\1/) {
			$hitDef = $2;
			$annot = ( $3 && $3 !~ /^\s*$/) ? $3 : "null";
		} elsif (/<(Hsp_bit-score)>(.*)<\/\1/) {
			$bits = $2;
		} elsif (/<(Hsp_evalue)>(.*)<\/\1/) {
			$evalue = $2;
		} elsif (/<(Hsp_query-from)>(.*)<\/\1/) {
			$qFrom = $2;
		} elsif (/<(Hsp_query-to)>(.*)<\/\1/) {
			$qTo = $2;
		} elsif (/<(Hsp_hit-from)>(.*)<\/\1/) {
			$hFrom = $2;
		} elsif (/<(Hsp_hit-to)>(.*)<\/\1/) {
			$hTo = $2;
		} elsif (/<(Hsp_query-frame)>(.*)<\/\1/){ # frame, added by Huang Fei, huangfei@genomics.cn
			$frame = $2;
		} elsif (/<(Hsp_identity)>(.*)<\/\1/) {
			$identity = $2;
		} elsif (/<(Hsp_gaps)>(.*)<\/\1/) {
			$gap = $2;
		} elsif (/<(Hsp_align-len)>(.*)<\/\1/) {
			$length = $2;
		} elsif (/<\/Hsp>/) {
			my $percent = sprintf("%.2f", $identity / $length * 100);
			my $hit = ($hitId =~ /gi/)? $hitId : $hitDef;
			my $s_e_format = ($frame > 0) ? "$qFrom\t$qTo" : "$qTo\t$qFrom";
			print OUT "$query\t$hit\t$percent\t$length\t" . ($length - $identity) . "\t$gap\t$s_e_format\t$hFrom\t$hTo\t$evalue\t$bits\t$annot\n";
			$gap = 0;
		}
	}
	close IN;

close OUT;

exit 0;
