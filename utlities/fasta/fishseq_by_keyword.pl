#!/usr/bin/perl
=head1 Name
	fishseq_by_keyword.pl

=head1 Description

	This perl is to extract the seq according to the seq_id's information,
	actually this is only a keyword in ID's, not keywords.

=head1 Contact & Version
		Author: Chentao YANG, yangchentao@genomics.org.cn
		Version: 1.0,  Date: 2014-12-1

=head1 Command-line Option
		perl   <infile | STDIN>
		--key <string>		the keyword you want
		--ex <string>	the keyword you do not want
		--help		output help information to screen

=head1 Usage Exmples
		perl  $0  -key ribosomal -ex protein >outfile

=cut
#################################################################################
use strict;
use warnings;
use Getopt::Long;
my ($Help,$Keyword,$except);
GetOptions(
	"help"=>\$Help,
	"key:s"=>\$Keyword,
	"ex:s"=>\$except
	);
#print "$Keyword";
die `pod2text $0` if (@ARGV == 0 ||$Help);
open IN,"$ARGV[0]" or die "Can not open file!";
$/=">";<IN>;$/="\n";
while (my $id=<IN>) {
	chomp $id;
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\s//g;
	$seq =~ s/\n//g;
	if ($Keyword && $except) {
		if ($id=~/$Keyword/ && !($id =~/$except/)){
			print ">$id\n" ;
			&fprint( $seq);
		}
	}
	if ($Keyword && !$except) {
		if ($id=~/$Keyword/) {
			print ">$id\n";
		#	&fprint ($seq);
			print "$seq\n";
		}
	}
	if ($except && !$Keyword){
		if (!($id =~/$except/)) {
			print ">$id\n" ;
			&fprint ($seq);
		}
	}
	$/="\n";
}
close IN;

sub fprint {
	my $sequence= shift;
	for ( my $pos = 0 ; $pos < length($sequence) ; $pos += 100 ) {
		print substr($sequence, $pos, 100)."\n";
	}
}
