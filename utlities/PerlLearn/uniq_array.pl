#!usr/bin/perl
my @array = ( 'a', 'a', 'a', 'a', 'a' );
my %count;
my @uniq_times = grep { ++$count{ $_ } < 2; } @array;
print "@uniq_times";
