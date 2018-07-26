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
		--t <output type> fa|fq output file type 
		--help		output help information to screen

=head1 Usage Exmples
		perl  ./simulate_NGS-reads.pl  -fa test.fa -l 150 -d 20 -i 300 

=cut

#################################################################################
use strict;
use warnings;
use Getopt::Long;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my ($n1,$n2) = ($mon+1,$mday);
#print "$n1,$n2\n";

my ($len,$fa,$depth,$insert,$type,$prefix,$Help);
GetOptions(
	"fa:s" => \$fa,
	"l:n" => \$len,
	"d:n" => \$depth,
	"i:n" => \$insert,
	"t:s" => \$type,
	help => \$Help
	);
die `pod2text $0` if ( !$fa || $Help);
$len ||= 150;
$depth ||= 30;
$insert ||= 300;
$type ||= "fq";
#$prefix ||= "K00133:316:HK3FJBBXX:$n1:$n2:";

$fa =~ /gz$/ ? open IN,"<:gzip",$fa : open IN ,"$fa" ;

my $file = $fa;
$file =~ s/.fa$//;
if ($type eq 'fq'){
	open L ,">$file"."_1.fq" or die "$!";
	open R ,">$file"."_2.fq" or die "$!";
}else {
	open L ,">$file"."_1.fa" or die "$!";
	open R ,">$file"."_2.fa" or die "$!";
}
my $quality = 'J'x $len;
my $count = 0;
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
					$count ++;
					$prefix = "K00133:316".":$id".":$n1:$n2".":$count:$count";
					if ($type eq 'fq'){
						print L "\@$prefix 1:N:0:NTGTTGGA\n$left\n+\n$quality\n";
						print R "\@$prefix 2:N:0:NTGTTGGA\n$right\n+\n$quality\n";
					
					}else {
						print L ">$id-$j/1\n$left\n";
						print R ">$id-$j/2\n$right\n";
					
					}
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


