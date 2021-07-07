#!/usr/bin/perl
=head1 Name
	find_com_diff.pl

=head1 Description
This perl is to find different or common lines between two files

=head1 Command-line Option
		perl   <infile | STDIN>
		--com
		--diff
		--help		output help information to screen

=head1 Usage Exmples
		perl  ./find_com_diff.pl  -diff -com file1(small) file2(large)

=cut

#################################################################################
use strict;
use warnings;
use Getopt::Long;
my ($diff,$com,$Help);
GetOptions(
	"diff" => \$diff,
	"com" => \$com,
	help => \$Help
	);
die `pod2text $0` if ( $Help || @ARGV ==0);
open SMA, shift;
open LAR, shift;
my %hash;
while (<SMA>) {
	chomp;
	$hash{$_} = 1;
}

if (defined $com){
	print "//Common\n";
	while (<LAR>) {
		chomp;
		print "$_\n" if ( exists $hash{$_});
	}
}

if (defined $diff) {
	print "//Diff\n";
	while (<LAR>) {
		chomp;
		print "$_\n" if ( ! exists $hash{$_});
	}
}
close SMA;
close LAR;
