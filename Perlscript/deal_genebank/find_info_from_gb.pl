#!/usr/bin/perl
use strict;
use warnings;
open F,"E:/sequence.gb";
my @a=<F>;
my $i=0;
my @accession;
my @mRNA_length;
my @description;
my $temp;
my $temp1;
my $temp2;
my $temp3;
my $genename;
my $locus_tag_name;
my @b;
my @CDS_length;
my @product;
my @protein_id;
my @AA_seq;
my @CDS_seq;
my @mRNA_seq;
my @gene_name;
my @tag_name;
my $start;
my $gene_id;
my @Gene_id;
while(i<=#a)
{
if(a[i]=~/^LOCUS(\s){7}/)
{
a[i]=~s/.*NM_(\d)+//;
a[i]=~s/(\s)+//g;
a[i]=~s/bp.*//s;
push @mRNA_length,a[i];$i++;
}elsif(a[i]=~/DEFINITION/)
{
a[i]=~s/DEFINITION  //;
a[i]=~s/\n//;
temp1=a[i];i++;
unless(a[i]=~/ACCESSION/)
{
a[i]=~s/(\s){11}//;
a[i]=~s/(\s)+$//;
a[i]=~s/\n//;
temp2=a[$i];
   temp1=temp1.$temp2;
   push @description,temp1;i++;
}else
{
push @description,temp1;i++;
}
temp1=;temp2='';
}
elsif(a[i]=~/VERSION/)
{
a[i]=~s/VERSION//;
a[i]=~s/GI.*//s;
       a[i]=~s/(\s)+//g;
push @accession,a[i];$i++;
}
elsif(a[i]=~/\/gene=/)
 {
 a[i]=~s/.*=//;
 a[i]=~s/\"//g;
 a[i]=~s/(\s)+//g;
 genename=a[i];i++;
     }elsif(a[i]=~/\/locus_tag=/)
{
        a[i]=~s/.*=//;
 a[i]=~s/\"//g;
 a[i]=~s/(\s)+//g;
 locustagname=a[i];i++;
}elsif(a[i]=~/\/db_xref=\"GeneID:/)
{
a[i]=~s/\/db_xref=\"GeneID://;
a[i]=~s/\"//g;
a[i]=~s/(\s)+//g;
geneid=a[i];i++;
}
elsif(a[i]=~/^(\s)+CDS(\s)+(\d)+/)
{
 a[i]=~s/CDS//;
         a[i]=~s/(\s)+//g;     
        a[i]=~s/\n//;
 @b=split/\.\./,a[i];
 temp=int(b[1])-int($b[0])+1;
 start=int(b[0]);
     push @CDS_length,$temp;
 @b=();$i++;
}
elsif(a[i]=~/\/product=\"/)
{
a[i]=~s/\/product=//;
a[i]=~s/\"//g;
a[i]=~s/^(\s)+//;
a[i]=~s/(\s)+$//;
a[i]=~s/\n//;
temp1=a[i];i++;
while(!(a[i]=~/\/protein_id=/))
{
a[i]=~s/^(\s){21}//;
a[i]=~s/(\s)+$//;
a[i]=~s/\n//;
temp2=a[$i];
  temp1=temp1.$temp2;
 $i++;
}
push @product,$temp1;
a[i]=~s/\/protein_id=//;
    a[i]=~s/\"//g;
   a[i]=~s/(\s)+//g;
   a[i]=~s/\n//;
push @protein_id,a[i];$i++;
   temp1=′′;temp2='';
}
elsif(a[i]=~/\/translation=/)
{
        a[i]=~s/\/translation=//;
a[i]=~s/\"//g;
a[i]=~s/^(\s)+//;
a[i]=~s/(\s)+$//;
a[i]=~s/\n//;
temp1=a[i];i++;
while(!(a[i]=~/ORIGIN/))
{
a[i]=~s/(\s)+//g;
a[i]=~s/\n//;
a[i]=~s/\"//;
temp2=a[$i];
  temp1=temp1.$temp2;
  $i++;
}
push @AA_seq,$temp1;     
    temp1=′′;temp2='';
}
elsif(a[i]=~/^ORIGIN(\s){6}/)
{
$i++;
   while(!(a[i]=~/\/\//))
{
a[i]=~s/(\s)+//g;
a[i]=~s/(\d)+//g;
a[i]=~s/\n//;
  temp2=a[$i];
  temp1=temp1.$temp2;
  $i++;
}
temp3=substr(temp1,start−1,temp);
push @CDS_seq,$temp3;
push @mRNA_seq,$temp1;
temp1=′′;temp2='';
}
elsif(a[i]=~/\/\//)
{
push @gene_name,$genename;
push @tag_name,$locus_tag_name;
push @Gene_id,$gene_id;
genename=′′;locus_tag_name='';$gene_id=''; 
$i++;
}
else
{
a[i]=~s/.*//sg;$i++;
}
}
open IN,">E:/seq_";
my $m;
my $str;
for(i=0;i<=#accession;i++)
{
print IN "Accession:\taccession[i]\nmRNA_length:\tmRNAlength[i]\n";
print IN "Description:\tdescription[i]\n";
print IN "Gene_id:\tGeneid[i]\nGene:\tgenename[i]\n";
print IN "Locus_tag:\t tagname[i]\nCDS_length:\tCDSlength[i]\n";
print IN "CDS_seq:\n";
$m=0;
$str='';
    while(m<=length(CDS_seq[$i]))
{
str=substr(CDS_seq[i],m,50);
m=m+50;
print IN "$str\n";
}
    print IN "Protein_id:\t proteinid[i]\n";
$m=0;
$str='';
print IN "AA_seq:\n";
while(m <= length(AA_seq[$i])) {
str=substr(AA_seq[i],m,50);
m=m+50;
print IN "$str\n";
}
$m=0;
$str='';
print IN "mRNA_seq:\n";
while(m<=length(mRNA_seq[$i]))
{
str=substr(mRNA_seq[i],m,50);
m=m+50;
print IN "$str\n";
}
print IN "\n";
}
