#!/usr/bin/perl -w
use strict;
use Bio::Perl;
print "Please input nucleotide sequence\n";
my $str = <STDIN> ;
chomp $str;
die "The sequence you inputed is not a standard cds" unless (((length $str) % 3) == 0);
my $translated = &TranslateDNASeq($str);

print "$translated\n";

sub TranslateDNASeq {
    use Bio::Seq;
	my $dna = shift;
	my $pep = translate_as_string($dna);
	return $pep;
}
