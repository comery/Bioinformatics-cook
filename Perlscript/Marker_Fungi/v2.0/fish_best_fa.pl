#!usr/bin/perl -w 
use strict;
die "perl $0 <*.fasta> <final.lst> <sample name>" unless (@ARGV == 3);
open IN,shift;#the file format is fasta;
open IN1,shift;#the file consist of the id of seq and fam_id;
my $tag = shift;
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
	my @yy=split /\s+/,$_;
	$bb{$yy[0]}=$yy[1];
}

foreach my $key (keys %bb) {
	print ">$key\_$bb{$key}_$tag\n$aa{$key}";
}

#while (my $fam_id=<IN2>) {
#	chomp $fam_id;
#	print ">$bb{$fam_id}_$fam_id\n$aa{$bb{$fam_id}}";
#}
