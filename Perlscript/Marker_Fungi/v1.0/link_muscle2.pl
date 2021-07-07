#!usr/bin/perl -w 
###################################################
#this perl is to change the format of muscle output
#both DNA and PRO seq can be accepted!
###################################################
use strict;
use List::Util qw(min max);
open IN,shift or die "$!!";
my ($xingxing,%hash,$id,$seq,%len,$id_len,$max_len_id);
while (<IN>){
	$_=~s/^\s+$//g;
	chomp;
	if ($_ eq ""||$_=~/^MUSCLE/|| $_=~/\*/) {
		next;
	}else{
		$id=(split /\s+/,$_)[0];
		$seq=(split /\s+/,$_)[1];
		$id_len=length ($id);
		$len{$id} = $id_len;
		if(!exists $hash{$id}){
			$hash{$id}=$seq;
		}else{
			$hash{$id}="$hash{$id}"."$seq";
		}
	}
}

################### print blank format ############
$max_len_id = max(values  %len);
my ($key,$value,@seq);
while (($key,$value)=each %hash){
	push @seq,$value;
	my $blank = " " x ($max_len_id - $len{$key} + 5);
	print "$key"."$blank"."$value\n";
}
print " " x ($max_len_id +5);

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
#my $xing_lens = &index($symble);
print "The whole length of aligned seq:\t$ll\nThe number of conserved bases:\t$conserved\nmax_repeat:\t$repeat\n";
my @blocks = &primer ($symble);
my $block_Num = @blocks;
print "potential blocks:\t$block_Num\n";
if ($block_Num >= 2 ) {
	my $left_block = min (@blocks);
	my $right_block = max (@blocks);
	my $exp_len = ($right_block - $left_block + 1) - 25;
	my $expected = substr($symble,$left_block,$right_block-$left_block+1);
	my $exp_con = ($expected =~ s/\*/\*/g);
	my $rate = $exp_con/$exp_len;
	print "The length of expected seq:\t$exp_len\n";
	print "The rate of conserved sites in expected seq:\t$rate";
}else {
	print "The length of expected seq:\tNA\n";
	print "The rate of conserved sites in expected seq:\tNA";
}
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
my @blocks;
my $n;
for ($n = 0;$n<=$len-25;$n++){
	my $temp_str = substr $str,$n,25;
	$count = ($temp_str =~ tr/*/*/);
	if ($count >=22){
		push @blocks,$n;
	}
}
return @blocks;
}
##############################
#sub index {
#	my $string = shift;
#	my $char = "*";
#	my $index_left = index($string,$char);
#	my $LL = reverse $string;
#	my $index_right = index($LL,$char);
#	my $xing_len =$ll-$index_right-$index_left+1;
#	return $xing_len;
#}
