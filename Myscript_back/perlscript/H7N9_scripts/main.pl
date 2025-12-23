#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my $usage = "Usage:\tperl $Bin/$0 <outbase.txt> <cluster_ref.fasta.cor.maskN.fa>\n";
die "$usage" unless (@ARGV == 2);

## STEP 1
open OUT1, ">step1.sh";
print OUT1 "perl $Bin/split_out.base.txt_based_on_ref_maskN.pl $file1 $file2 outdir";
`sh step1.sh`;

## TIP
print STDOUT "Tip:\n\tNow you have finished the first step, and next it need locate correctly,so it is necessary to use reference, the reference's GI number is bellow,download them with genebank format and save as reference.gb file. As all set, you can run the step2.sh\n";

if ( -f "id.info.list") {
	system ("cat id.info.list|awk '{print \$2}' ");
}else {
	die "some errors happened, please check the step1!";
}

## STEP 2
open OUT2, ">step2.sh";
print OUT2 <<LL;
if [ -f "reference.gb" ]; then
	perl $Bin/gb_cds_location.pl reference.gb >reference_cds_location.fa
else
	echo "Can not find file reference_cds_location.fa!"
fi

if [ -f "cds_location_in_segment.txt" ]; then
	rm cds_location_in_segment.txt
fi
	for((i=1;i<=8;i++))
	do 
		perl $Bin/locate.pl outdir/seg\$i/seg\$i-out.base.txt reference_cds_location.fa  >>cds_location_in_segment.txt
	done
	perl $Bin/makesh.pl cds_location_in_segment.txt >syn_nsyn.sh
	sh syn_nsyn.sh
LL

#`sh step2.sh`;
close OUT1;
close OUT2;
