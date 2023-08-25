#! /usr/bin/perl -w
use strict;
open IN,shift || die "Usage:\n\tperl $0 assembly.fa";
my ($scaffold,$numA,$numT,$numC,$numG,$numN,$sumA,$sumT,$sumC,$sumG,$sumN,$sumP,@scaffolds,$scaf,$len_scafs,$len_scaf);
$/=">";                       #����'>'�ָ��ļ�����ȡ��ÿ�ζ�ȡ������Ϊ����'>'֮�������
while(my $a=<IN>){        	      #���������ļ���������������$a
	chomp $a;      		       #�����з�ȥ��
	my @lines=split /\n/,$a;    #���ջ��з��ָ�
 	my $ID=shift @lines;        #ȡ������ >����һ�У��������е�ID
 	my $scaffold=join"",@lines;    #�ϲ������У���Ϊscaffold
    push @scaffolds,$scaffold;     #�����鱣��ÿ��scaffold
    $sumA +=($numA=$a=~s/A/A/g);    #ͳ��A�ĸ���
    $sumT +=($numT=$a=~s/T/T/g);
    $sumC +=($numC=$a=~s/C/C/g);
    $sumG +=($numG=$a=~s/G/G/g);
    $sumN +=($numN=$a=~s/N/N/g);
    $sumP=$sumA+$sumT+$sumC+$sumG+$sumN; #ͳ���ܸ���$sumP
}
$/="\n";
print"bases summary :\nA\t$sumA\nT\t$sumT\nC\t$sumC\nG\t$sumG\nN\t$sumN\n";  #�������

my @sortscaf= sort{length($b)<=>length($a)} @scaffolds;
#print @sortscaf;
$len_scafs=0;
while( $len_scafs<($sumP/2)){        #���ƴ���ܳ�scafsû�е�fasta�ļ��������ܳ�sumP��һ��
	$scaf=shift(@sortscaf);
	$len_scaf=length ($scaf);
	$len_scafs+=$len_scaf ;
 }
     print "The scaffold N50 is $len_scaf";     
close IN;
