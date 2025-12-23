#!usr/bin/perl -w
=head1 Description

	While you finished the annotation,you should get the unique sequence from scaffold sequence according to CDs file,
	This script is used to extract seq from CDS file ,,

=head1 Usage
		
	perl getMitoGscaf.pl <CDS> <fa> <outfilename>

=head1 Version
		
	Version :  1.0
	Created :  06/25/2014 09:28:12 PM

=head1 Contact

	Author       :	Shiyi Du
	Email        :	dushiyi@genomics.cn
	Organization :	BGI-NGB 

=cut

use strict;
use warnings;
die `pod2text $0` unless (@ARGV==3);
open CDS, $ARGV[0] or die"$ARGV[0] $!\n";
open FA, $ARGV[1] or die"$ARGV[1] $!\n";
open OUT, ">$ARGV[2]\_mitogeneScaf.fa" or die "$ARGV[2]\_mitogeneScaf.fa $!\n";
my %scaf;
$/ = ">";
<CDS>;
while(<CDS>) {
	chomp;
	my @a = split /\n/;
	my @b = split /locus=/, $a[0];
	my @c = split /:/, $b[1];
	$scaf{$c[0]} = 1;
}
close CDS;
<FA>;
while(<FA>) {
	chomp;
	my @A = split /\n/;
	my @B = split /\s+/, $A[0];
	if($scaf{$B[0]}) {
		print OUT ">$_";
	}
}
close FA;
close OUT;
