#!/usr/bin/perl -w 
use strict;
use warnings;
my $gbfile=shift;
open IN,"$gbfile";
open OUT,">GeneBank_ITS2_rRNA.fasta";
while (<IN>) {
    my ($E,$F);
    chomp;
    if(/^VERSION.+:(.+)/){    
        my $gi=$1;
        $/="ORIGIN";
        my $str = <IN>;
        $/="\n";
        my @line = split /\n/,$str;
        foreach my $i(0..$#line) {
            my $line = $line[$i];
            if ($line =~ /misc\_RNA.+?(\d+).+?(\d+)/) {
                my $start=$1;
                my $end=$2;
                my $line_next = $line[$i+1];
                if ($line_next=~ /internal\s+transcribed\s+spacer\s+2/) {
                    $F=$start;
                    $E=$end;
                }elsif (defined $line[$i+2] && $line[$i+2] =~ /internal\s+transcribed\s+spacer\s+2/) {
                    $F=$start;
                    $E=$end;
                }
            }
        }
        $/="//";
        my $list=<IN>;
        chomp $list;
        $/="\n";
        $list=~ s/\s+//g;
        $list=~ s/\d+//g;
        $list = uc($list);
        if (defined $F) {
            my $s28_len = $E-$F+1;
            my $fa=substr($list,$F-1,$s28_len);
            my $fa_len= length $fa;
            if ($fa_len >100) {
                print OUT ">gi:$gi\_ITS2\_$F\_$E\_$s28_len\n$fa\n";
            }
        }
    }
}
    close IN;
    close OUT;
