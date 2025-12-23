#!usr/bin/perl -w 
use strict;
use List::Util qw(min);
open IN,shift;
open OUT,">min_rate";
my $n=`awk 'BEGIN {FS=" "}{print NF}' test.ped|head -n 1 `;
my $num_bb=($n-6)/2;
while (<IN>) {
my @aa=split /\s+/,$_;
for (my $i=7;$i <=$n;$i+=2;my $j=1;$j <=$num_bb;$j++ ) {
my @$j;
push @$j,$aa[$i];
push @$j,$aa[$i+1];
}
for ($j=1;$j <=$num_bb;$j++){
	foreach (@$j){
	my $num_A=grep /A/g,@$j;
	my $num_T=grep /T/g,@$j;
	my $num_C=grep /C/g,@$j;
	my $num_G=grep /G/g,@$j;
	my $rate_A=$num_A/($num_A+$num_T+$num_C+$num_G);
	my $rate_T=$num_T/($num_A+$num_T+$num_C+$num_G);
	my $rate_C=$num_C/($num_A+$num_T+$num_C+$num_G);
	my $rate_G=$num_G/($num_A+$num_T+$num_C+$num_G);
	my @rate= ($rate_A,$rate_T,$rate_C,$rate_G); 
	my @new_rate;
		foreach $_(@rate) {
		push @new_rate,$_ if (!$_==0)
		}
		#print "@new_rate";
		my $min_rate = min(@new_rate);

		print "$rate_A\t$rate_T\t$rate_C\t$rate_G\t$min_rate\n";
		print OUT "$min_rate\n";
	}
}
close OUT;
`paste test.map min_rate >test.map.new`;

