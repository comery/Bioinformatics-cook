#!/usr/bin/perl

=head1 NAME

subsampleBysize.pl

=head1 SYNOPSIS

To get a specific size of data from a fastq file.

=head1 DESCRIPTION

Gzipped file is accepted.
The output data is extracted from the beginning of original file.
Accepts g, m, and k for easy input.

=head1 Contact

Author: Shangjin Tan, tanshangjin@genomics.cn

=head1 Version

Version: 1.0,  Date: 2019-01-11

=head1 OPTIONS

	--i	<string> 	the fastq file
	--s	<string> 	size to extract. Accpet g, m, and k. For example, 4g.
	--o	<string>	the output dir name,default dir is "./clean"
	--help	print this options

=head1 USAGE

perl subsampleBysize.pl -i Clean.fq.gz -s 4m -o out.gz

=head1 COPYRIGHT

Copyright 2019, Shangjin Tan.  All Rights Reserved.

=cut

##########################################################

use 5.010;
use warnings;
use strict;
use autodie;
use Getopt::Long;

use lib '/hwfssz1/ST_EARTH/Reference/ST_DIVERSITY/PUB/USER/tanshangjin/scripts/';
use MyModule;

my ($help, $in, $size, $out);
GetOptions(
	"i=s" => \$in,
	"s:s" => \$size,
	"o=s" => \$out,
	"help|h!" => \$help
);

die `pod2text $0` if ($help || !$size || !$in || !$out);
open OUT, "| gzip  > $out" or die $!;

if ($size =~ /g$/i){
	$size =~ s/g$//i;
	$size = $size * 1e9;
} elsif ($size =~ /m$/i){
	$size =~ s/m$//i;
	$size = $size * 1e6;
} elsif ($size =~ /k$/i){
	$size =~ s/k$//i;
	$size = $size * 1e3;
}
print "$size";
my $total = 0;
my $fh = open_file ($in);
while (<$fh>){
	my $id = $_;
	chomp (my $seq = <$fh>);
	$total += length $seq;
	my $third = <$fh>;
	my $qual = <$fh>;
	
	if ($total > $size){
		last;
	} else {
		print OUT "$id$seq\n$third$qual";
	}
}
close $fh;
close OUT;

exit;
