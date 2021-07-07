#!/usr/bin/perl
=head1 Name
According to the start and end sites to extract the substring

=head1 Command-line Option
	perl <infile|STDIN>
	--s <number>	the number of start site
	--e <number>	the number of end site
	--help		output this help information

=head1 Contact & Version
	Author: Chentao YANG, yangchentao@genomics.org.cn
	Version: 1.0,  Date: 2014-8-8

=head1 Usage Exmples
	perl substr.pl -s 10 -e 200 test.fa

=head1 Attention
	If start site number is larger than end site number, the strand is +,
	else, the strand is -; You should know that.

=cut

###############################################################
use strict;
use warnings;
use Getopt::Long;
my ($start,$end,$Help,);
GetOptions(
	"s:n"=>\$start,
	"e:n"=>\$end,
	"help:s"=>\$Help,
	)|| die "please use --help to get help\n";
die `pod2text $0` if ($Help || @ARGV == 0);

my $cut = $end - $start + 1;
print STDERR "cut length is $cut\n";

open IN,shift;
$/=">";<IN>;$/="\n";
my $seq_temp;
while(<IN>){
	my $id=$_;
	my $nid=(split /\s+/,$id)[0];
	chomp $nid;
    $/=">";
    my $seq=<IN>;
	$seq=~s/\n//g;
    $/="\n";
	my ($len,$nseq);
	if ($end > $start) {
		my $len=$end-$start+1;
		$nseq=substr($seq,$start-1,$len);
	}else{
		my $len=$start-$end+1;
		$seq_temp=substr($seq,$end-1,$len);
		$nseq=reverse($seq_temp);
	}
	print ">$nid:$start-$end\n";
	&warpper($nseq)
	#print "You are so careless with making the start equal the end \nso that you can just get only one site!" if ($len==1);
}

sub warpper {
	my $seq = shift;
	my $length=length $seq;
	for (my $i=0; $i<$length; $i+=60) {
		my $part=substr($seq,$i,60);
		print "$part\n";
	
	}
}


