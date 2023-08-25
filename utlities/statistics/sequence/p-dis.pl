#!usr/bin/perl -w 


my $param_string = <<__PARAMS__;

##---------------------------------------------------------------------#
This perl is to caculate the p-distance between two aligned sequence .
Usage: 
	
	perl $0 <test.tab>

the file test.tab is like :

	FARAN229-14.COI-5P      ---TCTTTACTTAATCTTTGGCGCTTGGGCCGGGATAGTAGGAACAGCCCTTAGCCTGCTCATTCGAGCAGAACTT
	EULEP1885-15.COI-5P     AACATTATATTTTATTTTTGGAATTTGAGCAGGAATAGTGGGAACATCTTTAAGAATCTTAATTCGAATAGAATTA

!!! THE SEQUENCE'S ID MUST BE UNIQUE!!!

__PARAMS__

use strict;
use List::Util qw(min max);
die "$param_string" unless (@ARGV == 1);
my $file = shift;
open IN,$file;
open OUT ,">$file.pdis.txt";
my ($xingxing,%hash);
while (<IN>){
	next if (/^MUSCLE/);
	next if (/^\s+$/);
	last if (/^>>symbols/);
	chomp;
	my @a = split;
	$hash{$a[0]} = uc($a[1]);
}

################### print blank format ############
my @head = keys %hash;
my (%distance,@distances);
while (@head) {
	my $que = shift (@head);
	foreach my $key (@head) {
	#	$distance{$que}{$key} = &dis($hash{$que},$hash{$key});
		my $diverg = &dis($hash{$que},$hash{$key});
		push @distances,$diverg;
		print "$diverg\n";
		print OUT "$que <----> $key\t$diverg\n";
	}
}
my $min_dis = min (@distances);
my $max_dis = max (@distances);
print  OUT "#The max p-distance is: $max_dis\n";
print  OUT "#The min p-distance is: $min_dis\n";
close IN;
close OUT;
#####################################################
sub dis {
	die "Pairwise model only accepts two values!\nUsage:\t\t my \$distance = &dis(\$str1,\$str2)\n" unless (@_ == 2);
	my @seq = @_;
	my $length = length ($_[0]);
	my ($y,@array);
	my $conserved = 0;
	my $gap = 0;
	my $vain = 0;
	for($y=0;$y<$length;$y++){
		foreach my $align(@seq){
			my $aa=substr $align,$y,1;
			push @array,$aa;
		}
		my %count;
		my @tmp=grep { ++$count{ $_ } < 2; } @array;
		if (@tmp==1 && $tmp[0]=~/[A-Z]/) {
	#		$symble.="*";
			$conserved++;
		}elsif(@tmp == 1 && $tmp[0] eq "-") { ##for removing the both gap(-) situation;
			$vain++;
		}elsif (@tmp ==2 && @tmp ~~ /-/){
			$gap++;
		}
		@array=();
	}
	print "$length\t$vain\t$gap\t$conserved\n";
	my $diverg = ($length-$vain-$gap-$conserved)/($length-$gap-$vain);
	return $diverg;
}
