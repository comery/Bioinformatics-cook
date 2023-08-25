#!/usr/bin/perl 
die "Usage: perl $0 <fasta>" unless (@ARGV == 1) ;

open IN,shift;

$/=">";<IN>;$/="\n";
while (my $id = <IN>) {
	chomp $id;
	$/= ">";
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/\n//g;
	print ">$id\n$seq\n";
	$/="\n";
	
}

close IN;
