#!/usr/bin/perl
use strict;
use Getopt::Long;
my ($list,$sub,$rev);
GetOptions(
	"list:s"=>\$list, #chr	start	end
	"sub:s"=>\$sub,#-sub chr:start:end
	"rev:s"=>\$rev,
);

my $usage = "Usages:
	perl $0 -list <list file> <fasta>
	perl $0 -sub <chr1:1:100> <fasta>
	perl $0 -rev <chr1:1:100> <fasta>";
if (@ARGV < 1){
	print $usage;
	exit()
}
my $fa=shift;
my %hash;
open IN,"$list";
while(<IN>){
	chomp;
	my @a=split /\s+/;
	$hash{$a[0]}{$a[1]}=$a[2];
}
close IN;

if(defined $sub){
	my ($temp_sa,$temp_st,$temp_en)=(split /:/,$sub)[0,1,2];
	$hash{$temp_sa}{$temp_st}=$temp_en;
}

if($fa=~/\.gz$/){
	open IN,"gzip -dc $fa |" || die "can't open the input file,$!";
}else{
	open IN,"$fa" || die "can't open the input file,$!";
}
$/=">";<IN>;
$/="\n";
while(<IN>){
	chomp;
	my $id=(split /\s+/)[0];
	$/=">";
	my $str=<IN>;
	$str =~ s/>$//g;
	$str =~ s/\n//g;
	my $str_len=length $str;
	if(exists $hash{$id}){
		foreach my $x(keys %{$hash{$id}}){
			my ($start,$end);
			if($x <= 0){$start =1;}else{$start=$x;}
			if($hash{$id}{$x} > $str_len){$end=$str_len;}else{$end=$hash{$id}{$x};}
			my $tem=substr ($str,$start-1,$end-$start+1);
			if(defined $rev){
				$tem = reverse $tem;
				$tem =~ tr/AGCTagct/TCGAtcga/;
			}
			print ">$id\_$start\_$end\n$tem\n";
		}
	}
	$/="\n";
}
close IN;

