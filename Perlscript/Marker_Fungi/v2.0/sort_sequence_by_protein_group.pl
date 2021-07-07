#!usr/bin/perl -w 
use strict;
die "Usage: perl $0 <*.fasta> <*.lst>" unless (@ARGV == 2);
open IN,shift;#the file format is fasta;
open IN1,shift;#the file consist of the id of seq and fam_id;
my (%aa,%bb,@yy,$fam_id,$id);
$/=">";<IN>;$/="\n";
while (<IN>) {
	chomp;
	$id=(split /\s+/,$_)[0];
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$/="\n";
	$aa{$id}=$seq;
}
while (<IN1>) {
	chomp;
	print ">$_\n$aa{$_}";
}

#while (my $fam_id=<IN2>) {
#	chomp $fam_id;
#	print ">$bb{$fam_id}_$fam_id\n$aa{$bb{$fam_id}}";
#}
