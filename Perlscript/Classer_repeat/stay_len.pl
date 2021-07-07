#!/usr/bin/perl -w
use strict;
my $file = shift;
open IN,"<$file" or die "$!";
my $sum ;
while  (<IN>) {
    chomp;
    my @a = split /\s+/,$_;
    my $len = $a[3] - $a[2];
 #   print "$len\n";
    $sum += $len;
  }
print "$file\t$sum\n";
