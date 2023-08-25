#!/usr/bin/perl -w

=pod
description: convert blast's result from -m0 (pair-wise) to -m8 (tabular, without comment lines)
author: Zhang Fangxian, zhangfx@genomics.cn
created date: 20090720
modified date: 20100410, 20090910, 20090909, 20090818, 20090727, 20090721
=cut

use Getopt::Long;

my ($input, $output, $help);

GetOptions(
	"input:s" => \$input,
	"output:s" => \$output,
	"help|?" => \$help,
);

sub usage {
	print STDERR << "USAGE";
description: convert blast's result from -m0 (pair-wise) to -m8 (tabular, without comment lines)
usage: perl $0 [options]
options:
	-input: input file, blast's result in format -m0
	-output: output tabular file, default is ./[-input].tab
	-help: help information
e.g.:
	perl $0 -input m0.txt -output m8.tab
USAGE
}

if (!defined $input || defined $help) {
	&usage();
	exit 1;
}

# check input file
if (!-f $input) {
	print "file $input not exists\n";
	exit 1;
}

# main
$output ||= &getFileName($input) . ".tab";

# iteration
open IN, "< $input" || die $!;
open OUT, "> $output" || die $!;

$flag = 0;
$noHitFlag = 0;
$itemIdx = 0;
$hitIdx = 0;
$hspIdx = 0;
while (<IN>) {
	chomp;
	if (/^BLAST/) {
		$itemIdx++;
		$hitIdx = 0;
		$flag = 1;
	} elsif (/^Reference:/) {
		$flag = 2;
	} elsif (/^Query=\s(.*)/) {
		$flag = 3;
		$_ = $1;
		$query = "";
	} elsif (/\s+\(([\d]+)\sletters\)$/) {
		$flag = 4;
		#$queryLen = $1;
	} elsif (/^Database:\s(.*)/) {
		$flag = 5;
		#$db = $1;
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
		#$length = $1;
		@temp = split(/ /, $hitDef, 2);
		$hitAccss = (split /\|/, $temp[0])[-1];
		$hitAccss =~ s/\..*$//;
	} elsif (/Score =\s+(.*)\sbits\s\(([\d]+)\),\sExpect.*=\s+(.*)/) {
		$hspIdx++;
		$flag = 8;
		$bits = $1;
		#$score = $2;
		$expect = $3;
	} elsif (/Identities =\s+([\d]+)\/([\d]+)\s.*,\s+Positives =\s+([\d]+)\/(.*)/) {
		$flag = 9;
		$identity = $1;
		$alignLen = $2;
		#$positive = $3;
		$gap = 0;
		if ($4 =~ /Gaps =\s+([\d]+)\//) {
			$gap = $1;
		}
	} elsif (/Frame =\s[+]*([\-\d]+)/) {
		$flag = 10;
		#$frame = $1;
	} elsif (/Query:\s(.*)/) {
		if ($flag != 11) {
			$flag = 11;
			@temp = split /\s+/, $1;
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
		#$reference .= "~" . $_ if ($_ ne "");
	} elsif ($flag == 3) {
		$query .= $_;
	} elsif ($flag == 6) {
		$hitDef .= " &gt;" . & trim($_);
		$hitDef = &trim($hitDef);
		$hitDef =~ s/^&gt;//;
	} elsif ($flag == 11) {
		if (/Query:\s(.*)/) {
			@temp = split /\s+/, $1;
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
			$noHitFlag = 0;
		}
	}

	if ($_ eq "") {
		if ($flag == 11 && $querySeq ne "") {
			print OUT (split /\s/, $query)[0] . "\t" . (split /\s/, $hitDef)[0] . "\t" . sprintf("%.2f", ($identity / $alignLen * 100)) . "\t$alignLen\t" . ($alignLen - $identity) . "\t$gap\t$queryFrom\t$queryTo\t$hitFrom\t$hitTo\t$expect\t$bits\n";
			$querySeq = "";
			$hitSeq = "";
			$midLine = "";
		}
	}
}

close OUT;
close IN;

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
