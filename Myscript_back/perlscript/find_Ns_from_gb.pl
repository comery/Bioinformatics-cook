#!/usr/bin/perl -w
use Bio::SeqIO;
my $inputfile = shift or die "$!";
my $in = Bio::SeqIO-> new(-file => "$inputfile", "-format" => 'genbank');
while (my $seq_obj=$in->next_seq()) {
	## get the whole complete mitochondrion/chloroplast genome sequence
	my $source = $seq_obj->seq;
	my $sou_len = $seq_obj->length;
	my $Ns;
	$source =~ /N/ ? $Ns = $source =~ s/N/N/g : $Ns = 0;
	my $len_withoutN = $sou_len - $Ns;
	print "$inputfile\t$sou_len\t$Ns\t$len_withoutN\n";
}