#!/usr/bin/perl -w 
use strict;
open IN, shift;
$/=">";<IN>;$/="\n";
while (my $id = <IN>) {
	chomp $id;
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	my $len = length ($seq);
	my ($i,$rev_complement);
	print ">$id";
	for ($i=0;$i<=$len;$i++){
		my $tem = substr($seq,$i,1);
		if ($tem eq "A"|| $tem eq "a") {
			$rev_complement .= "T";
		}elsif ($tem eq "T"||$tem eq "t"){
			$rev_complement .= "A";
		}elsif ($tem eq "G"||$tem eq "g"){
			$rev_complement .= "C" ;
		}elsif ($tem eq "C"||$tem eq "c") {
			$rev_complement .= "G";
		}else{
			$rev_complement .= "$tem" ;
		}
	}
	$rev_complement = reverse($rev_complement);
	print "$rev_complement\n";
	$/="\n";
}
