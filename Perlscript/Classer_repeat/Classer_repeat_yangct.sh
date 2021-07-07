#!/bin/bash
if [ $# -lt 1 ]
then 
     echo "sh $0 <gff>"
     exit
fi
gff=$1

perl -ne '@t=split/\t/;print "$t[0]\t"; print $1,"\t$t[3]\t$t[4]" if($t[8]=~/Class=([^;]+);/);print "\n"' $gff  | awk '($2~/^DNA/)'  >DNA.lst
perl -ne '@t=split/\t/;print "$t[0]\t"; print $1,"\t$t[3]\t$t[4]" if($t[8]=~/Class=([^;]+);/);print "\n"' $gff  | awk '($2~/^LINE/)' >LINE.lst
perl -ne '@t=split/\t/;print "$t[0]\t"; print $1,"\t$t[3]\t$t[4]" if($t[8]=~/Class=([^;]+);/);print "\n"' $gff  | awk '($2~/^SINE/)' >SINE.lst
perl -ne '@t=split/\t/;print "$t[0]\t"; print $1,"\t$t[3]\t$t[4]" if($t[8]=~/Class=([^;]+);/);print "\n"' $gff  | awk '($2~/^LTR/)'  >LTR.lst
/share/app/msort/bin/msort -k '{m1,m3}' DNA.lst    >DNA.lst.sort
/share/app/msort/bin/msort -k '{m1,m3}' LINE.lst   >LINE.lst.sort
/share/app/msort/bin/msort -k '{m1,m3}' LTR.lst    >LTR.lst.sort
/share/app/msort/bin/msort -k '{m1,m3}' SINE.lst   >SINE.lst.sort
perl  /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/perlscript/Classer_repeat/combine_repeat_by_type.pl DNA.lst.sort  DNA
perl  /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/perlscript/Classer_repeat/combine_repeat_by_type.pl LINE.lst.sort LINE
perl  /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/perlscript/Classer_repeat/combine_repeat_by_type.pl LTR.lst.sort  LTR
perl  /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/perlscript/Classer_repeat/combine_repeat_by_type.pl SINE.lst.sort SINE
perl -ne '@t=split /\t/;$hash{$t[1]} += $t[3]-$t[2]+1;if(eof) {foreach my $key (sort keys %hash) {print "$key\t$hash{$key}\n"}}' DNA.lst.sort  | msort -k '{m1,m3}' -  | sort -k2nr >DNA.lst.sort.class
perl -ne '@t=split /\t/;$hash{$t[1]} += $t[3]-$t[2]+1;if(eof) {foreach my $key (sort keys %hash) {print "$key\t$hash{$key}\n"}}' LINE.lst.sort | msort -k '{m1,m3}' -  | sort -k2nr >LINE.lst.sort.class
perl -ne '@t=split /\t/;$hash{$t[1]} += $t[3]-$t[2]+1;if(eof) {foreach my $key (sort keys %hash) {print "$key\t$hash{$key}\n"}}' LTR.lst.sort  | msort -k '{m1,m3}' -  | sort -k2nr >LTR.lst.sort.class
perl -ne '@t=split /\t/;$hash{$t[1]} += $t[3]-$t[2]+1;if(eof) {foreach my $key (sort keys %hash) {print "$key\t$hash{$key}\n"}}' SINE.lst.sort | msort -k '{m1,m3}' -  | sort -k2nr >SINE.lst.sort.class
perl /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/perlscript/Classer_repeat/stay_len.pl DNA.lst.sort  >lens.stat
perl /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/perlscript/Classer_repeat/stay_len.pl LINE.lst.sort >>lens.stat
perl /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/perlscript/Classer_repeat/stay_len.pl LTR.lst.sort  >>lens.stat
perl /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/perlscript/Classer_repeat/stay_len.pl SINE.lst.sort >>lens.stat
