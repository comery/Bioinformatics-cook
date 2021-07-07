#!usr/bin/perl -w
=head1 Name
	split_fastq.pl

=head1 Description
	this script is to split a big fastq file into two half files ;

=head1 Contact & Version
	Author: Chentao YANG, yangchentao@genomics.cn
	Version: 1.0,  Date: 2014-11-18

=head1 Command-line Option
	--fq1	*	<string> 	the *_1.gz of fastq file
	--fq2		<string>	the *_2.gz of fastq file [optional]
	--o		<string>	the output dir name,default dir is "./"
	--l		<number>	the length you want to trim from the end of reads.
	--help				print this options

=head1 Usage 
	perl split_fastq.pl -fq1 *_1.fq.gz -fq2 *_2.fq.gz -o  outdir -l 50

=cut
#######################################################
use strict;
use Getopt::Long;
use File::Basename;
#use PerlIO::gzip;
my ($FQ1,$FQ2,$Outdir,$Help,$trim);
GetOptions(
			"fq1:s" =>\$FQ1,
			"fq2:s" =>\$FQ2,
			"o:s" =>\$Outdir,
			"l:n" =>\$trim,
			"help" =>\$Help
			);
die `pod2text $0` if ($Help || ! $FQ1);
$Outdir ||= "./trimed";#default outdir is current dir
$trim ||= 50;
`mkdir  $Outdir` unless (-d $Outdir);
if ($FQ1 =~ /.gz$/) {
	open IN1,"gzip -dc $FQ1|" or die "Can't open $FQ1 !";
}else {
	open IN1,"$FQ1" or die "Can't open $FQ1 !";
}
if ($FQ2) {

	if ($FQ2 =~ /.gz$/) {
		open IN2,"gzip -dc $FQ2|" or die "Can't open $FQ2 !";
	}else {
		open IN2,"$FQ2" or die "Can't open $FQ2 !";
	}
	my $ffq2 = basename($FQ2);
	$ffq2 =~ s/.gz$//;
	open OUT2,"|gzip > $Outdir/$ffq2.gz";
	my $line2
}

#print "$len_l\t$len_r\n";
#die "The reads you give is too short to do triming!" unless ($len_l >= 100 && $len_r >= 100);

my $ffq1 = basename($FQ1);
$ffq1 =~ s/.gz$//;
open OUT1,"|gzip > $Outdir/$ffq1.gz";
my $line1;

while ($line1 = <IN1> ) {
	chomp $line1;
	my $id1 = &order($line1);
	my $seq1 = <IN1>;
	<IN1>;
	my $qual1 = <IN1>;
	chomp $seq1;chomp $qual1;
	substr($seq1,-$trim) = "";
	substr($qual1,-$trim) = "";
	
	print  OUT1 "\@$line1\n$seq1\n+\n$qual1\n";

	if ($FQ2){
		my $line2 = <IN2>;
		chomp $line2;
		my $id2 = &order($line2);
		die "fq1 and fq2 at are worry order! Check it out!" unless ($id1 eq $id2) ;
		my $seq2 = <IN2>;
		<IN2>;
		my $qual2 = <IN2>;
		chomp $seq2;chomp $qual2;
		substr($seq2,-$trim) = "";
		substr($qual2,-$trim) = "";
		print  OUT2 "\@$line2\n$seq2\n+\n$qual2\n";

	}

}
#	$/="\n";
close IN1;
close OUT1;
if ($FQ2) {
	close IN2;
	close OUT2;
}

#######################################################
###############  sub script ########################### 
#######################################################
sub order {
	my $str1 = shift;
#	my $lane1 = $1 if ($str1 =~ /(\S+)\/(\d)/);	
	my $lane1 = (split /\s+/,$str1)[0];
	return $lane1;
	}

