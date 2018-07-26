#!/usr/bin/perl -w 
use strict;
use warnings;
my $gbfile=shift;
open IN,"$gbfile";
open OUT,">GeneBank_5S_rRNA.fasta";
while (<IN>) {
    my ($E,$F,$gi);
    chomp;
    if(/^VERSION.+GI:(.+)/){    
        $gi=$1;
	}elsif(/^\s+ORGANISM\s+(.+)/){
		my $organism = $1;
        $/="ORIGIN";
        my $str = <IN>;
        $/="\n";
        my @line = split /\n/,$str;
        foreach my $i(0..$#line) {
            my $line = $line[$i];
            if ($line =~ /rRNA.+?(\d+).+?(\d+)/) {
                my $start=$1;
                my $end=$2;
                my $line_next = $line[$i+1];
                if ($line_next=~ /\b5S\b/) {
                    $F=$start;
                    $E=$end;
                }elsif (defined $line[$i+2] && $line[$i+2] =~ /\/.+\=\"\b5S\b/) {
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
            my $length = $E-$F+1;
            my $fa=substr($list,$F-1,$length);
           # my $fa_len= length $fa;
            if ($length >100) {
                print OUT ">gi|$gi| $organism 5S ribosomal RNA\n$fa\n";
            }
        }
    }
}
    close IN;
    close OUT;
