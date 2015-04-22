#!usr/bin/perl -w 
use strict;
open IN,shift;#.*blast.out
my (@a,$fam_id,$identity,$id1,$id2);
my $i=0;
while (<IN>){
	chomp;
	$identity=(split /\s+/,$_)[2];
	$id1=(split /\s+/,$_)[0];
	$id2=(split /\s+/,$_)[1];
	$fam_id=(split /_/,$id1)[-1];
	$i++ if ($identity>98);
}
print "$fam_id\n" if ($i<12);
close IN;
