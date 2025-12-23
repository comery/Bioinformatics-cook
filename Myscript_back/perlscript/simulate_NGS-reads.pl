#!/usr/bin/perl
=head1 Name
	simulate_NGS-reads.pl

=head1 Description
This perl is to simulate NGS reads by cutting fasta file randomly

=head1 Contact & Version
		Author: Chentao YANG, yangchentao@genomics.cn
		Version: 1.0,  Date: 2016-3-28

=head1 Command-line Option
		perl   <infile | STDIN>
		--fa <string>	fasta file
		--l <number>	read length you want cut into
		--d	<number>	average depth expected
		--i <number>	insert size length of your simulation
		--help		output help information to screen

=head1 Usage Exmples
		perl  ./simulate_NGS-reads.pl  -fa test.fa -l 150 -d 20 -i 300 

=cut

#################################################################################
use strict;
use warnings;
use Getopt::Long;
my ($len,$fa,$depth,$insert,$Help);
GetOptions(
	"fa:s" => \$fa,
	"l:n" => \$len,
	"d:n" => \$depth,
	"i:n" => \$insert,
	help => \$Help
	);
die `pod2text $0` if ( !$fa || $Help);
$len ||= 150;
$depth ||= 30;
$insert ||= 300;

$fa =~ /gz$/ ? open IN,"<:gzip",$fa : open IN ,"$fa" ;

my $file = $fa;
$file =~ s/.fa$//;

open L ,">$file"."_1.fa" or die "$!";
open R ,">$file"."_2.fa" or die "$!";

$/=">";<IN>;$/="\n";
while (<IN>) {
	chomp;
	my $id = $_;
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/ //g;
	$seq =~ s/\n//g;
	$seq= uc $seq;
	my @aa = &one2two ($seq);
		foreach my $i(@aa) {
			for my $j(1..$depth) {
				my $tmp = int(rand 20)-10+$insert;
				my $ll = $i;
				while (length $ll >= $tmp) {
					my $cut = substr $ll,0,$tmp;
					substr ($ll,0,$tmp) = "";
					my ($left,$right) = &pair($cut) ;
					print L ">$id-$j/1\n$left\n";
					print R ">$id-$j/2\n$right\n";
				}
			}
		}

	$/="\n"
}


sub one2two {
	my $str = shift;
	my $l = length $str;
	my $seed_A = int(rand $l);
	my $part1 = substr $str,0,$seed_A;
	substr ($str,0,$seed_A) = "";
	my @a = ($part1,$str);
	return @a;
}

sub pair {
	my $str = shift;
	my $left = substr $str,0,$len;
	my $revcom = reverse $str;
	$revcom =~ tr/ATCG/TAGC/;
	my $right = substr($revcom,0,$len);
	return $left,$right;
	}
