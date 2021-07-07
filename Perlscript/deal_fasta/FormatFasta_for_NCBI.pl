#!/usr/bin/perl -w
use strict;

die "Usage: perl $0 <fasta> <output file>" unless (@ARGV == 2) ;
my $in = $ARGV[0];
$in =~ /gz$/ ? open IN,"gzip -dc $in|" : open IN, "$in";

my $l = 200;
my $info = "[organism=Crocuta crocuta] [tech=wgs]";
my $line = 80;

print STDOUT "length cutoff is 200 bp\n";
print STDOUT "sequence information is \"$info\"\n";
print STDOUT "base in a line is 80 bp\n";
print STDOUT "processing...\n";

open OUT, ">$ARGV[1]";
$/=">";<IN>;$/="\n";
while (my $id = <IN>) {
	chomp $id;
	$id = $1 if ($id =~ /^(\S+)/);
	$/= ">";
	my $seq = <IN>;
	chomp $seq;
	my $len = &real_len($seq);
	if ($len >= $l) {
		print OUT ">$id $info\n";
		for (my $i=0; $i<$len; $i+=$line) {
			my $part=substr($seq,$i,$line);
			print OUT "$part\n";					
		}

	}
	$/="\n";
	
}

close IN;
close OUT;

sub real_len {
	my $str = shift;
	$str =~ s/[\s\n]//g;
	return length $str;
}

