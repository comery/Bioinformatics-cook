#!/usr/bin/perl
    $line ='ATCGGTACTGCTAGCTGCA';
  my @people=$line=~/\w{1}/g;
  my %count ;
  $count{$_}++ foreach @people;
  print %count;

