#!/usr/bin/perl
=head1 Name
	fishseq_by_keyword.pl

=head1 Description
this perl is to extract the seq according the seq_id's information,actually this is keyword in information,yes,this's keyword not keywords,if you have interesting in update it ,you can change this script!

=head1 Contact & Version
		Author: Chentao YANG, yangchentao@genomics.org.cn
		Version: 1.0,  Date: 2014-12-1

=head1 Command-line Option
		perl   <infile | STDIN>
		--key <string>		the keyword you want
		--ex <string>	the keyword you do not want
		--help		output help information to screen

=head1 Usage Exmples
		perl  ./fishseq_by_keyword.pl  -key ribosomal -ex protein >outfile

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
die `pod2text $0` if (@ARGV == 0 ||$Help);
open IN,shift;
$/=">";<IN>;$/="\n";
while (my $id=<IN>) {
	chomp $id;
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\s//g;
	if ($Keyword && $except) {
		print ">$id\n$seq\n" if ($id=~/$Keyword/ && !($id =~/$except/));
	}
	if ($Keyword && !$except) {
		print ">$id\n$seq\n" if ($id=~/$Keyword/);
	}
	if ($except && !$Keyword){
		print ">$id\n$seq\n" if (!($id =~/$except/));
	}
	$/="\n";
}
close IN;
