#! /usr/bin/perl -w
use strict;
open IN,shift || die "Usage:\n\tperl $0 assembly.fa";
my ($scaffold,$numA,$numT,$numC,$numG,$numN,$sumA,$sumT,$sumC,$sumG,$sumN,$sumP,@scaffolds,$scaf,$len_scafs,$len_scaf);
$/=">";                       #按照'>'分割文件并读取，每次读取的内容为两个'>'之间的内容
while(my $a=<IN>){        	      #按行输入文件，赋给标量变量$a
	chomp $a;      		       #将换行符去掉
	my @lines=split /\n/,$a;    #按照换行符分割
 	my $ID=shift @lines;        #取出含有 >的那一行，就是序列的ID
 	my $scaffold=join"",@lines;    #合并所有行，即为scaffold
    push @scaffolds,$scaffold;     #用数组保存每条scaffold
    $sumA +=($numA=$a=~s/A/A/g);    #统计A的个数
    $sumT +=($numT=$a=~s/T/T/g);
    $sumC +=($numC=$a=~s/C/C/g);
    $sumG +=($numG=$a=~s/G/G/g);
    $sumN +=($numN=$a=~s/N/N/g);
    $sumP=$sumA+$sumT+$sumC+$sumG+$sumN; #统计总个数$sumP
}
$/="\n";
print"bases summary :\nA\t$sumA\nT\t$sumT\nC\t$sumC\nG\t$sumG\nN\t$sumN\n";  #输出个数

my @sortscaf= sort{length($b)<=>length($a)} @scaffolds;
#print @sortscaf;
$len_scafs=0;
while( $len_scafs<($sumP/2)){        #如果拼接总长scafs没有到fasta文件中序列总长sumP的一半
	$scaf=shift(@sortscaf);
	$len_scaf=length ($scaf);
	$len_scafs+=$len_scaf ;
 }
     print "The scaffold N50 is $len_scaf";     
close IN;
