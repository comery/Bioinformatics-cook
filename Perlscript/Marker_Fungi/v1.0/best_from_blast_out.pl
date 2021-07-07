#!usr/bin/perl -w 
use  strict;
open IN,shift;
my (@a,@b,$fam_id);
while (my $line=<IN>) {
	chomp $line;
	@a=split /\s+/,$line;
	$fam_id=(split /_/,$a[1])[2];
	if (!grep /$a[0]/,@b) {
	print "$a[0]\t$fam_id\t$a[6]\t$a[7]\n";
	push @b,$a[0];
	}

}
	
