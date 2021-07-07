#!/usr/bin/perl
my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst) = localtime();    
$year += 1900;    
$mon++;
my $date = "$year-$mon-$day";    
print $date, "\n";
