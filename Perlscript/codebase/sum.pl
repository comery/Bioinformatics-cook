#!usr/bin/perl -w 
use strict;#count the number items whose value is zero!
my (@a,@b,$number,$sum);
open IN,shift;
open OUT,">sum";

while (<IN>) {
	@a=split /\s+/,$_;
	$number=grep /$a[0]/,@b;
	if ($number==0){
	$sum++;
	print OUT "$a[0]\n";
	push (@b,$a[0]);
	}
}
print $sum;
close IN;
close OUT;
