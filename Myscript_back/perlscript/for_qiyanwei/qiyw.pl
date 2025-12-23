#!usr/bin/perl -w 
use strict;
use List::Util qw(min);
open OUT,">min_rate";
my ($i,$n);
$n=`awk 'BEGIN {FS=" "}{print NF}' test.ped|head -n 1 `;
#print $n;
for ($i=7;$i <= $n;$i+=2) {

my $a=`awk '{print \$$i\$$i+1}' test.ped`;
$a=~s/\n//g;
#print $a;
my $num_A=($a=~s/A/A/g);
my $num_T=($a=~s/T/T/g);
my $num_C=($a=~s/C/C/g);
my $num_G=($a=~s/G/G/g);
#print "$num_T\n";
my $rate_A=$num_A/($num_A+$num_T+$num_C+$num_G);
my $rate_T=$num_T/($num_A+$num_T+$num_C+$num_G);
my $rate_C=$num_C/($num_A+$num_T+$num_C+$num_G);
my $rate_G=$num_G/($num_A+$num_T+$num_C+$num_G);
my @rate= ($rate_A,$rate_T,$rate_C,$rate_G); 
my @new_rate;
foreach $_(@rate) {
push @new_rate,$_ if (!$_==0)
}
print "@new_rate";
my $min_rate = min(@new_rate);

print "$rate_A\t$rate_T\t$rate_C\t$rate_G\t$min_rate\n";
print OUT "$min_rate\n";

}
close OUT;
`paste test.map min_rate >test.map.new`;

