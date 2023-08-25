#!usr/bin/perl -w
use strict;
use PerlIO::gzip;
my $file=shift;
open IN,"<:gzip","$file";
$/="@";<IN>;$/="\n";
while (my $a=<IN>) {
	my $b=<IN>;
	my $c=<IN>;
	$/="@";
	my $d=<IN>;
	chomp $d;
	my @array=split //,$d;
	my $len=length($d);
#	print "$len\n";
	my $n=0;
	foreach my $i(@array){
		my $y=ord($i)-64;
		if ($y<20) {
		$n++;
#		print "$n\n";
		}
	}
	print "\@$a$b$c$d" if (($n/$len)>0.1);
	$/="\n";
}
close IN;
