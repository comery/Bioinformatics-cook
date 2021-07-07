#!usr/bin/perl -w 
use strict;
open IN1,shift;#uniq.id.list
open IN,shift;
my (@a,$line,$fam_id,%hash);
while (<IN>) {
	@a=split;
	$line="$a[0]\t$a[3]\t$a[17]\t$a[18]";
	$fam_id=$a[3];
	$hash{$fam_id}=$line;
}

while (<IN1>) {
	chomp;
	print  "$hash{$_}\n";
}
close IN;
close IN1;

