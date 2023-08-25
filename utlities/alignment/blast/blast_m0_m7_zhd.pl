#!/usr/bin/perl -w

=pod
description: convert blast's result from -m0 (pair-wise) to -m7 (XML)
author: Zhang Fangxian, zhangfx@genomics.cn
created date: 20090720
modified date: 20100410, 20090818, 20090727, 20090721
modify by Deng Chaojin, dengzhj@genomics.cn, 2009-08-27, 12:30
=cut

use strict;
use Getopt::Long;

my ($input, $output, $group_volume, $cut_result, $num_suffix_len, $help);

GetOptions(
	"input:s"  => \$input,
	"output:s" => \$output,
	"gvol:i"   => \$group_volume,
	"mf"       => \$cut_result,
	"numlen:i" => \$num_suffix_len,
	"help|?"   => \$help,
);

sub usage {
	print STDERR << "USAGE";
description: convert blast's result from -m0 (pair-wise) to -m7 (XML)
usage: perl $0 [options]
options:
	-input <str>  input file, blast's result in format -m0
	-output <str> output xml file, default is ./[-input].xml
	-gvol <int>   group volume, that is, the number of <Iteration></Iteration> pair in one <BlastOutput><BlastOutput>
	-mf           cut the result *.xml into multi parts and output, one <BlastOutput><BlastOutput> in a file (*.xml.number)
	-numlen <int> the length of the number string suffix (default: 5)
	-help:        help information
e.g.:
	perl $0 -input m0.txt -output m7.xml
	perl $0 -gvol 1 -mf -numlen 5 -input gt_300_Px1.extend_scaf.fa.blast.nr -output test.xml
USAGE
}

if (!defined $input || defined $help) {
	&usage();
	exit 0;
}

# check input file
if (!-f $input) {
	print STDERR "file $input not exists\n";
	exit 1;
}

# main
$output ||= &getFileName($input) . ".xml";
$num_suffix_len ||= 5;

# param
my @temp = split /\n/, `tail -n 30 $input`;

my $index = 0;
for (0 .. $#temp) {
	$index = $_;
	last if ($temp[$index] =~ /^Matrix: /);
}

my $matrix = $temp[$index];
my $open = (split /,/, $temp[$index + 1])[0];
my $extend = (split /,/, $temp[$index + 1])[1];
my $expect = (split /:/, $temp[$index + 6])[0];

$matrix =~ s/Matrix:\s//;
$open =~ s/Gap Penalties: Existence:\s//;
$extend =~ s/\sExtension:\s//;
$expect =~ s/Number of sequences better than\s//;

# iteration
my ($flag, $noHitFlag, $hitDef, $t, $i, $line) = (0, 0, '', 0, 0, '');
my ($algorithm, $version, $reference, $db, $itemIdx, $query, $queryLen) = ('BLASTX', '2.2.18 [Mar-02-2008]', '~Reference:', 'nr ', 0, 'C5902753', 0);
my ($hitIdx, $hitAccss, $length) = (0, 'XP_002082351', 0);
my ($hspIdx, $bits, $score, $queryFrom, $queryTo, $hitFrom, $hitTo) = (1, 0, 0, 0, 0, 0, 0);
my ($frame, $identity, $positive, $gap, $alignLen, $querySeq, $hitSeq, $midLine) = (-1, 0, 0, 0, 0, 'MYRPN', 'MYRPN', 'HL   +  + SP');
my ($num_suffix_count, $num_suffix_str) = (1, '00001');
my $format_str = '%0' . $num_suffix_len . 'd';
my $out_part_file = "$output.$num_suffix_str";

open IN, "< $input" || die $!;
if ( !((defined $cut_result) and (defined $group_volume)) )
{ open OUT, "> $output" || die $!; }

if ( !((defined $cut_result) and (defined $group_volume)) )
{
	print OUT <<XML;
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
XML
}

while (<IN>) {
	chomp;
	if (/^BLAST/) {
		$itemIdx++;
		$hitIdx = 0;
		$flag = 1;
		($algorithm, $version) = split / /, $_, 2;
	} elsif (/^Reference:/) {
		$flag = 2;
		$reference = "";
	} elsif (/^Query=\s(.*)/) {
		$flag = 3;
		$_ = $1;
		$query = "";
	} elsif (/\s+\(([\d]+)\sletters\)$/) {
		$flag = 4;
		$queryLen = $1;
	} elsif (/^Database:\s*(\S*)/) {
		$flag = 5;
		$db = $1;

		if ( ($itemIdx-1) != 0 )
		{
			print OUT
"\t\t\t\t\t</Hit_hsps>\n\t\t\t\t</Hit>\n\t\t\t</Iteration_hits>\n\t\t</Iteration>\n";
		}

		if ( (defined $group_volume) and ( ($itemIdx-1) % $group_volume==0 ) and ( ($itemIdx-1) != 0 ) )
		{
			print OUT
"\t</BlastOutput_iterations>\n</BlastOutput>\n";
		}

		if (  ( (defined $group_volume) and ( ($itemIdx-1) % $group_volume==0 ) ) or ( !(defined $group_volume) and ( ($itemIdx-1)==0 ) )  )
		{
			if ( (defined $cut_result) and (defined $group_volume) )
			{
				if ( ($num_suffix_count-1)!=0 )
				{ close OUT; }

				$num_suffix_str = sprintf("$format_str", $num_suffix_count);
				$num_suffix_count++;
				$out_part_file = $output.'.'.$num_suffix_str;
				open OUT, "> $out_part_file" || die $!;

				print OUT <<XML;
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
XML
			}

			print OUT
'<BlastOutput>'."\n".
"\t<BlastOutput_program>$algorithm</BlastOutput_program>\n" .
"\t<BlastOutput_version>$algorithm $version</BlastOutput_version>\n" .
"\t<BlastOutput_reference>$reference</BlastOutput_reference>\n" .
"\t<BlastOutput_db>$db</BlastOutput_db>\n" .
"\t<BlastOutput_query-ID>$itemIdx</BlastOutput_query-ID>\n" .
"\t<BlastOutput_query-def>$query</BlastOutput_query-def>\n" .
"\t<BlastOutput_query-len>$queryLen</BlastOutput_query-len>\n" .
"\t<BlastOutput_param>\n" .
"\t\t<Parameters>\n" .
"\t\t\t<Parameters_matrix>$matrix</Parameters_matrix>\n" .
"\t\t\t<Parameters_expect>$expect</Parameters_expect>\n" .
"\t\t\t<Parameters_gap-open>$open</Parameters_gap-open>\n" .
"\t\t\t<Parameters_gap-extend>$extend</Parameters_gap-extend>\n" .
"\t\t</Parameters>\n" .
"\t</BlastOutput_param>\n" .
"\t<BlastOutput_iterations>\n";
		}

		print OUT
"\t\t<Iteration>\n" .
"\t\t\t<Iteration_iter-num>$itemIdx</Iteration_iter-num>\n" .
"\t\t\t<Iteration_query-ID>lcl|$itemIdx\_0</Iteration_query-ID>\n" .
"\t\t\t<Iteration_query-def>$query</Iteration_query-def>\n" .
"\t\t\t<Iteration_query-len>$queryLen</Iteration_query-len>\n" .
"\t\t\t<Iteration_hits>\n";
	} elsif (/^>(.*)/) {
		$hitIdx++;
		$hspIdx = 0;
		$flag = 6;
		$_ = $1;
		$hitDef = "";
		$querySeq = "";
		$hitSeq = "";
		$midLine = "";
	} elsif (/\s+Length =\s([\d]+)/) {
		$flag = 7;
		$length = $1;
		if ($hitIdx > 1) {
			print OUT <<XML;
\t\t\t\t\t</Hit_hsps>
\t\t\t\t</Hit>
XML
		}
		@temp = split(/ /, $hitDef, 2);
		$hitAccss = (split /\|/, $temp[0])[-1];
		$hitAccss =~ s/\..*$//;
		print OUT
"\t\t\t\t<Hit>\n" .
"\t\t\t\t\t<Hit_num>$hitIdx</Hit_num>\n" .
"\t\t\t\t\t<Hit_id>$temp[0]</Hit_id>\n" .
"\t\t\t\t\t<Hit_def>$temp[1]</Hit_def>\n" .
"\t\t\t\t\t<Hit_accession>$hitAccss</Hit_accession>\n" .
"\t\t\t\t\t<Hit_len>$length</Hit_len>\n" .
"\t\t\t\t\t<Hit_hsps>\n";

	} elsif (/Score =\s+(.*)\sbits\s\(([\d]+)\),\sExpect.*=\s+(.*)/) {
		$hspIdx++;
		$flag = 8;
		$bits = $1;
		$score = $2;
		$expect = $3;
	} elsif (/Identities =\s+([\d]+)\/([\d]+)\s.*,\s+Positives =\s+([\d]+)\/(.*)/) {
		$flag = 9;
		$identity = $1;
		$alignLen = $2;
		$positive = $3;
		$gap = 0;
		if ($4 =~ /Gaps =\s+([\d]+)\//) {
			$gap = $1;
		}
	} elsif (/Frame =\s[+]*([\-\d]+)/) {
		$flag = 10;
		$frame = $1;
	} elsif (/Query:\s(.*)/) {
		if ($flag != 11) {
			$flag = 11;
			@temp = split /\s+/, $1;
			if ($temp[0] > $temp[2]) {
				$t = $temp[0];
				$temp[0] = $temp[2];
				$temp[2] = $t;
			}
			$queryFrom = $temp[0];
			$queryTo = $temp[2];
			$querySeq = $temp[1];
			$i = index($_, $temp[1]);
			chomp($line = <IN>);
			$midLine = substr($line, $i);
			chomp($line = <IN>);
			@temp = split /\s+/, $line;
			$hitFrom = $temp[1];
			$hitTo = $temp[3];
			$hitSeq = $temp[2];
			<IN>;
			next;
		}
	} elsif (/.*No hits found/) {
		$noHitFlag = 1;
		$flag = 12;
	} elsif (/\s+Database:/) {
		$flag = 13;
	}

	if ($flag == 2) {
		s/\"/&quot;/;
		$reference .= "~" . $_ if ($_ ne "");
	} elsif ($flag == 3) {
		$query .= $_;
	} elsif ($flag == 6) {
		$hitDef .= " &gt;" . & trim($_);
		$hitDef = &trim($hitDef);
		$hitDef =~ s/^&gt;//;
	} elsif ($flag == 11) {
		if (/Query:\s(.*)/) {
			@temp = split /\s+/, $1;
			if ($temp[0] > $temp[2]) {
				$t = $temp[0];
				$temp[0] = $temp[2];
				$temp[2] = $t;
			}
			$queryTo = $temp[2];
			$querySeq .= $temp[1];
			$i = index($_, $temp[1]);
			chomp($line = <IN>);
			$midLine .= substr($line, $i);
			chomp($line = <IN>);
			@temp = split /\s+/, $line;
			$hitTo = $temp[3];
			$hitSeq .= $temp[2];
			<IN>;
		}
	} elsif ($flag == 12) {
		if ($noHitFlag == 1) {
			print OUT
"\t\t\t\t<Hit>\n" .
"\t\t\t\t\t<Hit_hsps>\n";
			$noHitFlag = 0;
		}
	}

	if ($_ eq "") {
		if ($flag == 11 && $querySeq ne "") {
			print OUT
"\t\t\t\t\t\t<Hsp>\n" .
"\t\t\t\t\t\t\t<Hsp_num>$hspIdx</Hsp_num>\n" .
"\t\t\t\t\t\t\t<Hsp_bit-score>$bits</Hsp_bit-score>\n" .
"\t\t\t\t\t\t\t<Hsp_score>$score</Hsp_score>\n" .
"\t\t\t\t\t\t\t<Hsp_evalue>$expect</Hsp_evalue>\n" .
"\t\t\t\t\t\t\t<Hsp_query-from>$queryFrom</Hsp_query-from>\n" .
"\t\t\t\t\t\t\t<Hsp_query-to>$queryTo</Hsp_query-to>\n" .
"\t\t\t\t\t\t\t<Hsp_hit-from>$hitFrom</Hsp_hit-from>\n" .
"\t\t\t\t\t\t\t<Hsp_hit-to>$hitTo</Hsp_hit-to>\n" .
"\t\t\t\t\t\t\t<Hsp_query-frame>$frame</Hsp_query-frame>\n" .
"\t\t\t\t\t\t\t<Hsp_identity>$identity</Hsp_identity>\n" .
"\t\t\t\t\t\t\t<Hsp_positive>$positive</Hsp_positive>\n" ;
			print OUT "\t\t\t\t\t\t\t<Hsp_gaps>" . $gap . "</Hsp_gaps>\n" if ($gap);
			print OUT
"\t\t\t\t\t\t\t<Hsp_align-len>$alignLen</Hsp_align-len>\n" .
"\t\t\t\t\t\t\t<Hsp_qseq>$querySeq</Hsp_qseq>\n" .
"\t\t\t\t\t\t\t<Hsp_hseq>$hitSeq</Hsp_hseq>\n" .
"\t\t\t\t\t\t\t<Hsp_midline>$midLine</Hsp_midline>\n" .
"\t\t\t\t\t\t</Hsp>\n";
			$querySeq = "";
			$hitSeq = "";
			$midLine = "";
		}
	}
}
print OUT "\t\t\t\t\t</Hit_hsps>\n\t\t\t\t</Hit>\n\t\t\t</Iteration_hits>\n\t\t</Iteration>\n\t</BlastOutput_iterations>\n</BlastOutput>\n";
close OUT;

exit 0;

sub trim {
	my ($s) = @_;
	$s =~ s/^\s+//;
	$s =~ s/\s+$//;
	return $s;
}

sub getFileName {
	my ($file_name) = @_;
	$file_name = (split /[\/\\]/, $file_name)[-1];
	$file_name =~ s/\.[^\.]*$//;
	return $file_name;
}
