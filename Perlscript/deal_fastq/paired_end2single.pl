#!usr/bin/perl -w
=head1 Name
	split_fastq.pl

=head1 Description
	this script is to merge paired-end fastq file into one file by order ;

=head1 Contact & Version
	Author: Chentao YANG, yangchentao@genomics.cn
	Version: 1.0,  Date: 2014-11-18

=head1 Command-line Option
	--fq1		<string> 	the *_1.gz of fastq file
	--fq2		<string>	the *_2.gz of fastq file 
	--o			<string>	the output dir name,default dir is "./"
	--help				print this options

=head1 Usage 
	perl $0 -fq1 *_1.fq.gz -fq2 *_2.fq.gz -o  outdir

=cut
#######################################################
use strict;
use Getopt::Long;
use File::Basename;
use PerlIO::gzip;
my ($FQ1,$FQ2,$Outfile,$Help);
GetOptions(
			"fq1:s" =>\$FQ1,
			"fq2:s" =>\$FQ2,
			"o:s" =>\$Outfile,
			"help" =>\$Help
			);
die `pod2text $0` if ($Help || ! $FQ1);

open IN1,"<:gzip", "$FQ1" or die "Can't open $FQ1 !";
open IN2,"<:gzip", "$FQ2" or die "Can't open $FQ2 !";
open OUT,">:gzip", "$Outfile";
my ($id1,$id2,$qual1,$qual2,$seq1,$seq2);
my $count=0;

while ($id1 = <IN1>) {
	$id2 = <IN2>;
	$seq1 = <IN1>;
	$seq2 = <IN2>;
	<IN1>;
	<IN2>;
	$qual1 = <IN1>;
	$qual2 = <IN2>;
	print  OUT "$id1$seq1+\n$qual1$id2$seq2+\n$qual2";
}
close IN1;
close IN2;


