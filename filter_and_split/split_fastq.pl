#!usr/bin/perl -w
=head1 Name
	split_fastq.pl

=head1 Description
	this script is to split a big fastq file into two half files ;

=head1 Contact & Version
	Author: Chentao YANG, yangchentao@genomics.cn
	Version: 1.0,  Date: 2014-11-18

=head1 Command-line Option
	--fq1		<string> 	the *_1.gz of fastq file
	--fq2		<string>	the *_2.gz of fastq file 
	--o			<string>	the output dir name,default dir is "./"
	--help				print this options

=head1 Usage 
	perl split_fastq.pl -fq1 *_1.fq.gz -fq2 *_2.fq.gz -o  outdir

=cut
#######################################################
use strict;
use Getopt::Long;
use File::Basename;
use PerlIO::gzip;
my ($FQ1,$FQ2,$Outdir,$Help);
GetOptions(
			"fq1:s" =>\$FQ1,
			"fq2:s" =>\$FQ2,
			"o:s" =>\$Outdir,
			"help" =>\$Help
			);
die `pod2text $0` if ($Help || ! $FQ1);
$Outdir ||= "./";#default outdir is current dir
`mkdir -p $Outdir/split_1 $Outdir/split_2`;
open IN1,"<:gzip", "$FQ1" or die "Can't open $FQ1 !";
open IN2,"<:gzip", "$FQ2" or die "Can't open $FQ2 !";
my $ffq1 = basename($FQ1);
my $ffq2 = basename($FQ2);
open OUT11,">:gzip", "$Outdir/split_1/$ffq1";
open OUT12,">:gzip", "$Outdir/split_1/$ffq2";
open OUT21,">:gzip", "$Outdir/split_2/$ffq1";
open OUT22,">:gzip", "$Outdir/split_2/$ffq2";
my ($line1,$line2);
my $count=0;
$/="@";<IN1>;<IN2>;
while ($line1 = <IN1>) {
	$line2 = <IN2>;
	chomp $line1;
	chomp $line2;
	my $id1 = &order($line1);
	my $id2 = &order($line2);
#	print "$id1\n$id2\n";
	if ($id1 ne $id2) {
		die "fq1 and fq2 at are worry order! Check it out!";
	}
	$count++;
	if ($count % 2 == 0) {
		print  OUT21 "\@$line1";
		print  OUT22 "\@$line2";
	}else{
		print  OUT11 "\@$line1";
		print  OUT12 "\@$line2";
	}
}
$\="\n";
close IN1;
close IN2;
close OUT11;
close OUT12;
close OUT21;
close OUT22;
#######################################################
###############  sub script ########################### 
#######################################################
sub order {
	my $str1 = shift;
	my $lane1 = $1 if ($str1 =~ /(\S+)\/(\d)/);
	return $lane1;
	}
	










