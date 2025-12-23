use strict;    
use warnings;    
    
    
my $dna      ='';    
my $protein  ='';    
my @file_data=( );    
my @filedata;  
my $revcom='';  
  
  
#打开文件   
@filedata  = get_file_data();  
#得到序列  
$dna       = extract_sequence_from_fasta_data(@filedata);    
  
#六框阅读翻译  
  
print "\n---------------------Reading Frame 1-----------------\n";  
$protein=translate_frame($dna,1);  
print_sequence($protein,70);  
  
print "\n---------------------Reading Frame 2-----------------\n";  
$protein=translate_frame($dna,2);  
print_sequence($protein,70);  
  
print "\n---------------------Reading Frame 3-----------------\n";  
$protein=translate_frame($dna,3);  
print_sequence($protein,70);  
  
print "\n---------------------Reading Frame 4-----------------\n";  
$protein=translate_frame($dna,4);  
print_sequence($protein,70);  
  
print "\n---------------------Reading Frame 5-----------------\n";  
$protein=translate_frame($dna,5);  
print_sequence($protein,70);  
  
print "\n---------------------Reading Frame 6-----------------\n";  
$protein=translate_frame($dna,6);  
print_sequence($protein,70);  
  
sub get_file_data  
{    
    # A subroutine to get data from a file given its filename  
    #读取文件的子序列  
    my $dna_filename;  
    my @filedata;  
    print "please input the Path just like this f:\\\\perl\\\\data.txt\n";     
    chomp($dna_filename=<STDIN>);   
    open(DNAFILENAME,$dna_filename)||die("can not open the file!");      
    @filedata     = <DNAFILENAME>;    
    close DNAFILENAME;    
    return @filedata;#子函数的返回值一定要记住写  
}  
  
sub extract_sequence_from_fasta_data    
{    
    #*******************************************************************    
    # A subroutine to extract FASTA sequence data from an array    
    # 得到其中的序列    
    # fasta格式介绍：    
    # 包括三个部分    
    # 1.第一行中以>开头的注释行，后面是名称和序列的来源    
    # 2.标准单字母符号的序列    
    # 3.*表示结尾    
    #*******************************************************************    
    
    my (@fasta_file_data) =@_;    
    my $sequence =' ';    
    foreach my $line (@fasta_file_data)    
    {    
        #这里忽略空白行    
        if ($line=~/^\s*$/)    
        {    
            next;    
        }    
        #忽略注释行    
        elsif($line=~/^\s*#/)    
        {    
            next;    
        }    
        #忽略fasta的第一行    
        elsif($line=~/^>/)    
        {    
            next;    
        }    
        else    
        {    
            $sequence .=$line;    
        }    
    }    
    $sequence=~s/\s//g;    
    return $sequence;    
}    
    
sub print_sequence    
{    
    # A subroutine to format and print sequence data    
    my ($sequence, $length) = @_;    
    for (my $pos =0; $pos<length($sequence);$pos+=$length)    
    {    
        print substr($sequence,$pos,$length),"\n";    
    }    
}    
    
       
    
sub codon2aa       
{       
    
    #第三种方法      
    #也就是运用哈希      
    #我们将所有的密码子作为hash的key，然后将代表的氨基酸作为hash的value      
    #然后进行匹配      
    # codon2aa       
    # A subroutine to translate a DNA 3-character codon to an amino acid       
    # Version 3, using hash lookup       
    my($codon) = @_;       
       
    $codon = uc $codon;#uc=uppercase;lc=lowercase      
                   #也就是大小写转换，uc表示将所有的小写 转换为大写      
               #lc将所有的大写转换为小写      
        
    my(%genetic_code) = (       
           
    'TCA' => 'S',    # Serine       
    'TCC' => 'S',    # Serine       
    'TCG' => 'S',    # Serine       
    'TCT' => 'S',    # Serine       
    'TTC' => 'F',    # Phenylalanine       
    'TTT' => 'F',    # Phenylalanine       
    'TTA' => 'L',    # Leucine       
    'TTG' => 'L',    # Leucine       
    'TAC' => 'Y',    # Tyrosine        
    'TAT' => 'Y',    # Tyrosine       
    'TAA' => '_',    # Stop       
    'TAG' => '_',    # Stop       
    'TGC' => 'C',    # Cysteine       
    'TGT' => 'C',    # Cysteine       
    'TGA' => '_',    # Stop       
    'TGG' => 'W',    # Tryptophan       
    'CTA' => 'L',    # Leucine       
    'CTC' => 'L',    # Leucine       
    'CTG' => 'L',    # Leucine       
    'CTT' => 'L',    # Leucine       
    'CCA' => 'P',    # Proline       
    'CCC' => 'P',    # Proline       
    'CCG' => 'P',    # Proline       
    'CCT' => 'P',    # Proline       
    'CAC' => 'H',    # Histidine       
    'CAT' => 'H',    # Histidine       
    'CAA' => 'Q',    # Glutamine       
    'CAG' => 'Q',    # Glutamine       
    'CGA' => 'R',    # Arginine       
    'CGC' => 'R',    # Arginine       
    'CGG' => 'R',    # Arginine       
    'CGT' => 'R',    # Arginine       
    'ATA' => 'I',    # Isoleucine       
    'ATC' => 'I',    # Isoleucine       
    'ATT' => 'I',    # Isoleucine       
    'ATG' => 'M',    # Methionine       
    'ACA' => 'T',    # Threonine       
    'ACC' => 'T',    # Threonine       
    'ACG' => 'T',    # Threonine       
    'ACT' => 'T',    # Threonine       
    'AAC' => 'N',    # Asparagine       
    'AAT' => 'N',    # Asparagine       
    'AAA' => 'K',    # Lysine       
    'AAG' => 'K',    # Lysine       
    'AGC' => 'S',    # Serine       
    'AGT' => 'S',    # Serine       
    'AGA' => 'R',    # Arginine       
    'AGG' => 'R',    # Arginine       
    'GTA' => 'V',    # Valine       
    'GTC' => 'V',    # Valine       
    'GTG' => 'V',    # Valine       
    'GTT' => 'V',    # Valine       
    'GCA' => 'A',    # Alanine       
    'GCC' => 'A',    # Alanine       
    'GCG' => 'A',    # Alanine       
    'GCT' => 'A',    # Alanine           
    'GAC' => 'D',    # Aspartic Acid       
    'GAT' => 'D',    # Aspartic Acid       
    'GAA' => 'E',    # Glutamic Acid       
    'GAG' => 'E',    # Glutamic Acid       
    'GGA' => 'G',    # Glycine       
    'GGC' => 'G',    # Glycine       
    'GGG' => 'G',    # Glycine       
    'GGT' => 'G',    # Glycine       
    );       
       
    if(exists $genetic_code{$codon})       
    {       
        return $genetic_code{$codon};       
    }      
    else      
    {       
       
            print STDERR "Bad codon \"$codon\"!!\n";       
            exit;       
    }       
}       
    
sub dna2peptide    
{    
    my ($dna)=@_;    
    my $protein ='';    
    for (my $i=0; $i<(length($dna)-2);$i+=3)    
    {    
        $protein .=codon2aa(substr($dna,$i,3));    
    }    
    return $protein;
}    
  
sub translate_frame  
{  
    my ($seq,$start,$end)=@_;  
    my $protein;  
      
    unless($end)  
    {  
        $end=length($seq);  
    }  
    return dna2peptide(substr($seq,$start-1,$end-$start+1));  
}  