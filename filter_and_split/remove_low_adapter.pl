#!usr/bin/perl -w
=head1 Name
	remove_low_adapter.pl

=head1 Description
	this script is to filter low quality reads and containing adapter reads.

=head1 Contact & Version
	Author: Chentao YANG, yangchentao@genomics.cn
	Version: 1.0,  Date: 2014-11-18

=head1 Command-line Option
	--fq1		<string> 	the *_1.gz of fastq file
	--fq2		<string>	the *_2.gz of fastq file 
	--o		<string>	the output dir name,default dir is "./clean"
	--a 	<str>	1.adapter.list.gz
	--b 	<str>	2.adapter.list.gz
	--q 	QUAL,CONT    filter out low quality reads, 'B'<->66
	--help				print this options

=head1 Usage
	perl remove_low_adapter.pl -fq1 *_1.fq.gz -fq2 *_2.fq.gz -a 1.adapter.list.gz -b 2.adapter.list.gz -o  outdir ./cleand -q 10,50 

=cut
#######################################################
use strict;
use Getopt::Long;
use File::Basename;
#use PerlIO::gzip;
my ($FQ1,$FQ2,$Outdir,$ad1,$ad2,$qual,$N,$Help,$trim,);
my $adapter = 0;
my $low = 0;
my $nn = 0;
my $total = 0;
my $clean = 0;
GetOptions(
			"fq1=s" =>\$FQ1,
			"fq2=s" =>\$FQ2,
			"o:s" =>\$Outdir,
			"a=s" =>\$ad1,
			"b=s" =>\$ad2,
			"q:s" =>\$qual,
			"n:s" =>\$N,
			"help" =>\$Help
			);
die `pod2text $0` if ($Help || ! $FQ1 || !$FQ2 || !$ad1 || !$ad2);
$Outdir ||= "./clean";#default outdir is current dir
$qual ||= "10,50";
$N ||= 10;
my $zl = (split /,/,$qual)[0];
my $cont = (split /,/,$qual)[1];

`mkdir  $Outdir` unless (-d $Outdir);
open LOG,">$Outdir/log";

#open adapter list
if ($ad1 =~ /.gz$/) {
	open IN3,"gzip -dc $ad1|" or die "Can't open $ad1 !";
}else {
	open IN3,"$ad1" or die "Can't open $ad1 !";
}
if ($ad2 =~ /.gz$/) {
	open IN4,"gzip -dc $ad2|" or die "Can't open $ad2 !";
}else {
	open IN4,"$ad2" or die "Can't open $ad2 !";
}

my (%ada1,%ada2);
while (<IN3>) {
	chomp;
	next if (/^#reads_id/);
	my @a = split /\s+/,$_;
	$ada1{$a[0]} =1 if ($a[3] < 90);
}

while (<IN4>) {
	chomp;
	next if (/^#reads_id/);
	my @a = split /\s+/,$_;
	$ada2{$a[0]} =1 if ($a[3] < 90);
}

close IN3;
close IN4;

#open fastq file
if ($FQ1 =~ /.gz$/) {
	open IN1,"gzip -dc $FQ1|" or die "Can't open $FQ1 !";
}else {
	open IN1,"$FQ1" or die "Can't open $FQ1 !";
}
if ($FQ2 =~ /.gz$/) {
	open IN2,"gzip -dc $FQ2|" or die "Can't open $FQ2 !";
}else {
	open IN2,"$FQ2" or die "Can't open $FQ2 !";
}

my $ffq1 = basename($FQ1);
my $ffq2 = basename($FQ2);
$ffq1 =~ s/.gz$//;
$ffq2 =~ s/.gz$//;
open OUT1,"|gzip > $Outdir/clean_$ffq1.gz";
open OUT2,"|gzip > $Outdir/clean_$ffq2.gz";

my ($line1,$line2);
my $rm_ada1 = 0;
my $rm_ada2 = 0;
$/="@";<IN1>;<IN2>;$/="\n";
while ($line1 = <IN1> ) {
	$line2 = <IN2>;
	$total++;
	chomp $line1;chomp $line2;
	my $id1 = &order($line1);
	my $id2 = &order($line2);
#	print "$id1\t$id2\n";
	if ($id1 ne $id2) {
		die "fq1 and fq2 are at worry order! Check it out!";
	}
	$/="@";

	#read lines left:
	my $str1 = <IN1>;my $str2 = <IN2>;
	chomp $str1;chomp $str2;
	#print "$str2\n";
	my @aa = split /\n/,$str1;
	my @bb = split /\n/,$str2;
	my $seq1 = $aa[0];my $qual1 = $aa[2];
	my $seq2 = $bb[0];my $qual2 = $bb[2];
	
	if (exists $ada1{$id1} ){
		$rm_ada1++;
		delete $ada1{$id1};
		last if (! %ada1); # if the %ada1 is empty so that stop to check!
	}elsif (exists $ada2{$id2}) {
		$rm_ada2++;
		delete $ada2{$id2};
		last if (! %ada2); # if the %ada2 is empty so that stop to check!
	}elsif (! &cut_q($qual1) || ! &cut_q($qual2)) {
		$low ++;	#count the low quality reads number (exactlly pairs);
	}elsif ( &cut_n($seq1)>$N || &cut_n($seq2)>$N ) {
		$nn ++;
	}else{
		$clean ++;
		print  OUT1 "\@$line1\n$seq1\n+\n$qual1\n";
		print  OUT2 "\@$line2\n$seq2\n+\n$qual2\n";
	}
	
	$/="\n";

}
#statistic info
my $is_adapter = $rm_ada1 + $rm_ada2;
print LOG "total reads pairs:\t$total\n";
print LOG "is_adapter pairs:\t$is_adapter\n";
print LOG "low_$qual:\t$low\n";
print LOG "N_$N:\t$nn\n";
print LOG "clean:\t$clean\n";

close IN1;
close IN2;
close OUT1;
close OUT2;
close LOG;

#######################################################
###############  sub script ########################### 
#######################################################
sub order {
	my $str1 = shift;
#	my $lane1 = $1 if ($str1 =~ /(\S+)\/(\d)/);	
	my $lane1 = (split /\s+/,$str1)[0];
	return $lane1;
	}

sub cut_n {
	my $str = shift;
	my $n = $str =~ tr/N/N/;
	return $n;
}

sub cut_q {
	my $str = shift;
	my $len = length $str;
	#print "$len\n";
	my $bad = 0;
	for (my $i=0;$i<$len;$i++){
		my $tmp = substr($str,$i,1);
		my $ord = ord($tmp);
		$bad++ if ($ord <$zl);
	}
	my $bad_rate = $bad/$len;
	$bad_rate < ($cont/100) ? return 1 :return 0;
}

