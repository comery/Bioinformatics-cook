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
	print ">$key\n$value\n";
}
