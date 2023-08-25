#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;
use Getopt::Long;
=head1 Name
statistic the information of assembly result(fasta)

=head1 Contact & Version
  Author: Chentao YANG, yangchentao@genomics.cn
  Version: 1.0,  Date: 2016-4-29
  
=head1 Usage Exmples
  perl  ./$0 assembly.fasta -l
  
=cut

#------------------------------------------------------------#
my ($out_len,$help);
GetOptions(
	"l" => \$out_len,
	"help" => \$help
	) || die "Please use --help option to get help\n";

die `pod2text $0` if ($help );
die "usage: perl $0 [ contig file ] [-l]\n" unless (@ARGV >= 1);

my (%hash,$base,$base_freeN,$aa,$gg,$cc,$tt,$nn,$frags);
open IN,$ARGV[0] or die "$!\n";

if ($ARGV[0] =~ /.gz$/) {
	open IN," gzip -dc $ARGV[0] | ";
}else{
	open IN,$ARGV[0];
}

my($head,$seq);
my @line;
$/=">";<IN>;$/="\n";
my (@length, @seq_length);
my $l_100 = 0 ;
my $l_500 =0 ;
my $l_1k = 0;
my $l_10k = 0;
my $l_100k = 0;
my $l_1M = 0 ;
open OUT,">$ARGV[0].length" if ($out_len);
while (<IN>) {
	chomp;
	@line=split;
    my $head=$line[0];        
    $/=">";
    my $seq=<IN>;
    chomp $seq;
	$seq=~s/[\s\n]//g;
   	my $len=length($seq);
  	if ($len >= 1000000) {
  		$l_1M ++;
  	}elsif($len >= 100000){
  		$l_100k ++;
  	}elsif($len >= 10000){
  		$l_10k ++;
  	}elsif($len >= 1000){
  		$l_1k ++;
  	}elsif($len >= 500){
  		$l_500 ++;
  	}elsif($len >= 100){
  		$l_100 ++;
  	}
	my $frag = &break($seq);
	$frags += $frag;
   	$base += $len;
	push @length,$len;
	my $AA = ($seq =~ s/A/A/gi);
	my $CC = ($seq =~ s/C/C/gi);
	my $GG = ($seq =~ s/G/G/gi);
	my $TT = ($seq =~ s/T/T/gi);
	my $NN = ($seq =~ s/N/N/gi);
	$aa += $AA;
	$cc += $CC;
	$gg += $GG;
	$tt += $TT;
	$nn += $NN;
	my $seq_len = $len - $NN;
	$base_freeN += $seq_len;
	
	print OUT "$head\t$len\t$seq_len\n" if ($out_len);
    $/="\n";
}

#GC_Content
my $GC_Content = ($cc + $gg)/$base; 
#longest,shortest,mean,median size
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@length);
my $shortest = $stat->quantile(0); # => 1
my $lower_quartile = $stat->quantile(1); # => 3.25	
my $median = $stat->quantile(2); # => 5.5
my $upper_quartile = $stat->quantile(3); # => 7.75
my $longest = $stat->quantile(4); # => 10
my $Scaffold_Num = @length;
my $mean = $base/$Scaffold_Num;

#about N

my $Average_length_of_N_in_scaffold = $nn/$Scaffold_Num;
#N50
my @sortscaf= sort{$b<=>$a} @length;
#make a sorted length hash
my $hash;
for my $m (0..$#sortscaf){
	$hash{$m} = $sortscaf[$m];
}
my $scaf = 0;

my ($N10,$N20,$N30,$N40,$N50,$N60,$N70,$N80,$N90);
for my $n (0..$#sortscaf){
#	my $len_scaf = shift @sortscaf;
	$scaf += $hash{$n};
	$N10 = $hash{$n+1} if ($scaf < 0.1*$base);
	$N20 = $hash{$n+1} if ($scaf < 0.2*$base);
	$N30 = $hash{$n+1} if ($scaf < 0.3*$base);
	$N40 = $hash{$n+1} if ($scaf < 0.4*$base);
	$N50 = $hash{$n+1} if ($scaf < 0.5*$base);
	$N60 = $hash{$n+1} if ($scaf < 0.6*$base);
	$N70 = $hash{$n+1} if ($scaf < 0.7*$base);
	$N80 = $hash{$n+1} if ($scaf < 0.8*$base);
	$N90 = $hash{$n+1} if ($scaf < 0.9*$base);

}


#output
open STA,">$ARGV[0].Statistics" or die "$!\n";
print STA "<-- Information for assembly Scaffold '$ARGV[0]'.(cut_off_length < 100bp) -->\n";
print STA "Size_includeN\t$base\n";
print STA "Size_withoutN\t$base_freeN\n";
print STA "Scaffold_Num\t$Scaffold_Num\n";
print STA "Mean_Size\t$mean\n";
print STA "Median_Size\t$median\n";
print STA "Longest_Seq\t$longest\n";
print STA "Shortest_Seq\t$shortest\n";
print STA "Singleton_Num";
print STA "Average_length_of_break(N)_in_scaffold\t$Average_length_of_N_in_scaffold\n";
print STA "scaffolds>100\t$l_100\n";
print STA "scaffolds>500\t$l_500\n";
print STA "scaffolds>1K\t$l_1k\n";
print STA "scaffolds>10K\t$l_10k\n";
print STA "scaffolds>100K\t$l_100k\n";
print STA "scaffolds>1M\t$l_1M\n";
print STA "Nucleotide_A\t$aa\n";
print STA "Nucleotide_C\t$cc\n";
print STA "Nucleotide_G\t$gg\n";
print STA "Nucleotide_T\t$tt\n";
print STA "GapN\t$nn\n";
printf STA "GC_Content\t$GC_Content\t(G+C)/(A+C+G+T)\n";

printf STA "N10\t$N10\n";
printf STA "N20\t$N20\n";
printf STA "N30\t$N30\n";
printf STA "N40\t$N40\n";
printf STA "N50\t$N50\n";
printf STA "N60\t$N60\n";
printf STA "N70\t$N70\n";
printf STA "N80\t$N80\n";
printf STA "N90\t$N90\n";

close OUT;
close STA;
#---------sub----------------------#
sub break {
	my $str = shift;
	if ($str =~ /N/) {
		my @a = split /N+/,$str;
		my $frags = @a;
		return $frags;
	}else {
		return 1;
	}
}

