#!usr/bin/perl -w
my $Usage = qq(
	perl remove_low_adapter.pl *.fastq(*.fq.gz) [number]
	
);
#######################################################
use strict;
use File::Basename;
die "$Usage" unless (@ARGV==2);
my $FQ = $ARGV[0];
my $remainder = $ARGV[1];
die "your number is not a int!\n"  unless ($remainder =~ m/^\d+$/);

#open fastq file
if ($FQ =~ /.gz$/) {
	open IN,"gzip -dc $FQ|" or die "Can't open $FQ !";
}else {
	open IN,"$FQ" or die "Can't open $FQ!";
}
my $ffq = basename($FQ);
$ffq =~ s/.gz$//;
open OUT,"|gzip > adrink_$ffq.gz";

my $total = 0;
my $keep =0;
my $line;
while ($line = <IN> ) {
	chomp $line;
	$total ++;
	my $seq = <IN>;
	chomp $seq;
	my $seq_len = length $seq;
	<IN>;
	my $Qual = <IN>;
	chomp $Qual;
	if ($total%$remainder == 0){

		print OUT "$line\n$seq\n+\n$Qual\n" ;
		$keep ++;
	}
}

print "total reads number:\t$total\n";
print "keep reads number:\t$keep\n";



