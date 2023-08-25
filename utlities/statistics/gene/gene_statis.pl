#!/usr/bin/perl -w

=head1 Name 

	gene_statis.pl

=head1 Description

	The file which was read by gene_statics.pl contains data of gene information of annotation. 
	The programm will calculate the information of gene, such as the Average length of Exon.

=head1 Version

	Authors: Zijun Xiong, xiongzijun@genomics.cn
	Version: 1.0, Date: 2010-04-18

=head1 Usuage

	perl gene_statis.pl <inputfile1> <inputfile2> 
	<inputfile1>: the format of file must be gff 
	<inputfile2>: the file must be the file which contains the information of cds

=head1 Example 

	perl gene_statis.pl ./../gene2.complete.gff ./../gene2.complete.gff.cds 

=cut

use strict;

die `pod2text $0` if (scalar (@ARGV) != 2);

my $inputfile1 = shift @ARGV; 
my $inputfile2 = shift @ARGV;

##calculate the Average length of cds and the other information###
open F,"<$inputfile1" or die "can't open file $0:$!";

my ($mrna_Count,$Total_CDSLen,$Aver_CDSLen,$last,$Intron_Len,$Intron_num);
my $CDS_Len=0;
my ($mrna_Len,$len,$Aver_IntronLen,$Aver_ExonCount,$Aver_ExonLen);
my $Aver_rnaLen;
my $CDS_Count;
while(<F>) {
	if(/mRNA\s+(\w+)\s+(\w+)/) {
		$len=$2-$1+1;
		$mrna_Len += $len;
		$mrna_Count++; 
		$Total_CDSLen += $CDS_Len;
		$last = 0;
		$CDS_Len = 0;
	}
	if(/CDS\s+(\d+)\s+(\d+)/) {
		$CDS_Count++;
		$CDS_Len += $2-$1+1;	
		if ($last and $last+1<$1){
			$Intron_Len += $1-$last-1;
			$Intron_num++;
		}
		$last = $2;
	}
}
close F;

$Total_CDSLen += $CDS_Len;
$Aver_CDSLen=$Total_CDSLen/$mrna_Count;
$Aver_IntronLen=$Intron_Len/$Intron_num;
$Aver_rnaLen=$mrna_Len/$mrna_Count;
$Aver_ExonCount=$CDS_Count/$mrna_Count;
$Aver_ExonLen=$Total_CDSLen/$CDS_Count;

###calcutate the Average GC content###
open IN,"<$inputfile2" or die "$!";
my $gc_count;
my $Aver_GC = 0;
$/=">";
<IN>;
$/="\n";
while(<IN>) {
	chomp(my $name = $_);
	$/=">";
	my $seq = <IN>;
	$/="\n";
	$seq =~ s/>$//;
	$seq =~ s/[\r\n]//g;
	my $gc = 0;
	for my $i (0..length($seq)-1) {
		$gc++ if(substr($seq,$i,1)) =~ /[GC]/i;
	}
	$gc_count += $gc/length($seq);
}
close IN;
$Aver_GC = $gc_count/$mrna_Count;


print "Total genes\tAver gene len(bp)\tAver CDS len(bp)\tAver CDS GC Ratio\tAver Exons per gene\tAver Exon len(bp)\tAver Intron len(bp)\n";
print "$mrna_Count\t$Aver_rnaLen\t$Aver_CDSLen\t$Aver_GC\t$Aver_ExonCount\t$Aver_ExonLen\t$Aver_IntronLen\n";
exit 0;
