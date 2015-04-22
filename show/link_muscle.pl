#!usr/bin/perl -w 
###################################################
#this perl is to change the format of muscle output
#both DNA and PRO seq can be accepted!
###################################################
use strict;
use List::Util qw(max);
open IN,shift or die "$!!";
my ($xingxing,%hash,$id,$seq,$id_len);
while (<IN>){
	$_=~s/^\s+$//g;
	chomp;
	if ($_ eq ""||$_=~/^MUSCLE/|| $_=~/\*/) {
		next;
	}else{
		$id=(split /\s+/,$_)[0];
		$seq=(split /\s+/,$_)[1];
		$id_len=length ($id);
		if(!exists $hash{$id}){
			$hash{$id}=$seq;
		}else{
			$hash{$id}="$hash{$id}"."$seq";
		}
	}
}

my ($key,$value,@seq);
while (($key,$value)=each %hash){
	push @seq,$value;
	print "$key\t$value\n";
}
print " " x $id_len ."\t";

my $length=length($seq[0]);
my ($y,@array,$symble);
my $conserved=0;
for($y=0;$y<$length;$y++){
	foreach my $align(@seq){
		my $aa=substr $align,$y,1;
		push @array,$aa;
	}
	my %count;
	my @tmp=grep { ++$count{ $_ } < 2; } @array;
	if (@tmp==1 && $tmp[0]=~/[A-Z]/) {
		$symble.="*";
		$conserved++;
	}else{
		$symble.="-";
	}
	@array=();
}
print "$symble\n";
my $repeat=&repeat_region($symble);
my $ll = length ($symble);
my $out = &primer ($symble);
my $xing_lens = &index($symble);
my $rate = $conserved/$xing_lens;
print "length\t\t$ll\nConserved\t$conserved";
close IN;

###############sub##############
sub repeat_region{
my $lens=shift;
my $long=length ($lens);
my  $j=0;my @count;my $i;
for ($i=0;$i<=$long-1;$i++){
	my $base=substr $lens,$i,1;
	if ($base eq "*"){
	$j++;
	}else{
	push @count,$j;
	$j=0;
	}
}
my $max_repeat=max(@count);
return $max_repeat;
}	
###############################
sub primer {
my $str = shift;
my $len = length($str);
my $count;
my $blocks = 0;
my $n;
for ($n = 0;$n<=$len-21;$n+=20){
	my $temp_str = substr $str,$n,20;
	$count = ($temp_str=~tr/*/*/);
	if ($count >=12){
		$blocks++;
	}
}
return $blocks;
}
##############################
sub index {
	my $string = shift;
	my $char = "*";
	my $index_left = index($string,$char);
	my $LL = reverse $string;
	my $index_right = index($LL,$char);
	my $xing_len =$ll-$index_right-$index_left+1;
	return $xing_len;
}
