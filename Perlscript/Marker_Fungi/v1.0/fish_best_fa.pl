#!usr/bin/perl -w 
use strict;
open IN,shift;#the file format is fasta;
open IN1,shift;#the file consist of the id of seq and fam_id;
open IN2,"com_all_best_to_best.id";
my (%aa,%bb,@yy,$fam_id,$id);
$/=">";<IN>;$/="\n";
while (<IN>) {
	chomp;
	$id=$_;
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$/="\n";
	$aa{$id}=$seq;
}
while (<IN1>) {
	chomp;
	my @yy=split /\s+/,$_;
	$bb{$yy[1]}=$yy[0];
}
while (my $fam_id=<IN2>) {
	chomp $fam_id;
	print ">$bb{$fam_id}_$fam_id\n$aa{$bb{$fam_id}}";
}
