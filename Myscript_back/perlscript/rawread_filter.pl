#! usr/bin/perl -w
use strict;
use Getopt::Long;

=head1 Description

        Usage: perl readfilter.pl <parameter>

        -lis    raw data path including fq and adpter, in the format of fq1\nfq2\nadpter1\nadapter2\n
        -par    how many task run in parallel, 2 by default
	-nma 	the maximun N allowed in single end reads, by default eq 10
	-qua 	the maixmun low quality (B) allowed in single end reads, by default eq 20
        -mis 	the cutoff of mismatch number of adapter, by default eq 3
	-aln 	the cutoff the length of alignment, by default eq 15
	-zip 	the output file compressed or not, (y | n), by defalt eq y;
	-help   print out this information

=cut

my ($Lis,$Par,$Help,$Nma,$Qua,$Mis,$Aln,$Zip);
GetOptions(
        "lis:s"=>\$Lis,
        "par:i"=>\$Par,
	"nma:i"=>\$Nma,
	"qua:i"=>\$Qua,
	"mis:i"=>\$Mis,
	"aln:i"=>\$Aln,
	"zip:s"=>\$Zip,
        "help"=>\$Help
);
die `pod2text $0` if ($Help || !defined ($Lis));

open LO, ">rawdatafilter.log" || die $!;

my $line=`less -S $Lis |wc -l`;
die "the list should contain both end of fastq and adapter\n"
unless ($line%4 == 0);

$Par ||=2;
$Nma ||=10;
$Qua ||=20;
$Mis ||=3;
$Aln ||=15;
$Zip = "y" if (!defined $Zip);
open IN, $Lis or die $!;
my @array;
while(<IN>){
        chomp;
        chomp(my $seq2=<IN>);
        chomp(my $adapter1=<IN>);
        chomp(my $adapter2=<IN>);
        push @array, [$_, $seq2, $adapter1, $adapter2];
}
close IN;

for (my $ii=0;$ii<=$Par-1;$ii++){
        my $pid = fork();
        if ($pid < 0){
                print "errors happened during pid generation\n";
        }elsif($pid == 0){
                filter ($ii);
                print "$ii has been finished";
                exit;
        }else{
                warn "the $ii process is waiting on line\n";
        }
}

sub filter {
	my $id =shift (@_);
	for (my $i=$id;$i<=$#array;$i+=$Par){
		my @a=split /\//,$array[$i][0];
                my @b=split /\//,$array[$i][1];
                my @c=split /\//,$array[$i][2];
                my @d=split /\//,$array[$i][3];
		die "wrong order" 
                unless ($a[-2] eq $b[-2] && $b[-2] eq $c[-2] && $d[-2] eq $c[-2]);
                `mkdir $a[-2]`;
		chdir "./$a[-2]/";
		my ($adf,$quf,$nf)=proc ($array[$i][0],$array[$i][1],$array[$i][2],$array[$i][3],$a[-2]);
		chdir "../";
		print LO "$a[-2]\nadapter contamination removed $adf\nlow quality removed $quf\nN nmuber removed $nf\n";
	}
}
close LO;
sub proc {
	my %hash;
	my $fq1=shift(@_);
	my $fq2=shift(@_);
	my $ad1=shift(@_);
	my $ad2=shift(@_);
	my $out=shift(@_);
	my $out1="$out\.1.fq";
	my $out2="$out\.2.fq";
	if ($fq1=~/gz$/){open (FFQ, "<:gzip",$fq1) }else{open (FFQ, $fq1)}
	if ($fq2=~/gz$/){open (RFQ, "<:gzip",$fq2) }else{open (RFQ, $fq2)}
	if ($ad1=~/gz$/){open (FAD, "<:gzip",$ad1) }else{open (FAD, $ad1)}
	if ($fq2=~/gz$/){open (RAD, "<:gzip",$ad2) }else{open (RAD, $ad2)}
	if ($Zip eq "y"){
		open (OUF, ">:gzip","$out1.gz");
		open (OUR, ">:gzip","$out2.gz")
	}else{
		open OUF, ">$out1";
		open OUR, ">$out2"
	}
	FIR:while(<FAD>){
		next if (/^#/);
		my @a=split;
		next FIR if ($a[-1]>$Mis and $a[-2]<$Aln);
		my $id=(split /\#/,$a[0])[0];
		$hash{$id}=1;
	}
	close FAD;
	SEC:while(<RAD>){
		next if (/^#/);
		my @a=split;
		next SEC if ($a[-1]>$Mis and $a[-2]<$Aln);
		my $id=(split /\#/,$a[0])[0];
		$hash{$id}=1;
	}
	close RAD;
	my $adf=0;
	my $quf=0;
	my $nf=0;
	THI:while(my $tf=<FFQ>){
		chomp($tf);
		chomp(my $tr=<RFQ>);
		chomp(my $seqf=<FFQ>);
		chomp(my $seqr=<RFQ>);
		chomp(my $orif=<FFQ>);
		chomp(my $orir=<RFQ>);
		chomp(my $quaf=<FFQ>);
		chomp(my $quar=<RFQ>);
		my $id=(split /\#/,$tf)[0];
		$id=~s/\@//;
		if (exists $hash{$id}){
			$adf++;
			next THI	
		}
		my $nseqf=$seqf=~s/N/N/g;
		my $nseqr=$seqr=~s/N/N/g;
		if ($nseqf>=$Nma or $nseqr>=$Nma){
			$nf++;
			next THI
		}
		my $nquaf=$quaf=~s/B/B/g;
		my $nquar=$quar=~s/B/B/g;
		if ($nquaf>=$Qua or $nquar>=$Qua){
			$quf++;
			next THI
		}
		print OUF "$tf\n$seqf\n$orif\n$quaf\n";
		print OUR "$tr\n$seqr\n$orir\n$quar\n";
	}
	close FFQ;
	close RFQ;
	close OUF;
	close OUR;
	return ($adf,$quf,$nf);
}
