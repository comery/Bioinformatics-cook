#!usr/bin/perl -w 
###################################################
#this perl is to change the format of muscle output
#both DNA and PRO seq can be accepted!
###################################################
use strict;
use List::Util qw(min max);
open IN,shift or die "$!!";
my (%hash,$id,$seq,%len,$id_len,$max_len_id);
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
#1####################################################
my $length=length($seq[0]);
my ($y,@array);
my $consensus ;
for($y=0;$y<=$length-1;$y++){
	foreach my $align(@seq){
		my $aa=substr $align,$y,1;
		push @array,$aa;
	}
	my %count;
	my @tmp=grep { ++$count{ $_ } < 2; } @array;
	if (@tmp==1 && $tmp[0]=~/[ATCG]/) {
		$consensus.=$tmp[0];
	}elsif (grep /-/,@tmp) {
		$consensus.="-";
	}elsif (@tmp==2 && $tmp[0]=~/[ATCG]/ && $tmp[1]=~/[ATCG]/){
		my $spot = &hybird(@tmp);
		$consensus.=$spot;
	}elsif (@tmp>2) {
		$consensus.="-";
	}
	@array=();
}
print "$consensus\n";

#2###############################################################
my %location;
my $n;
my $len = length ($consensus);
my $tag =0;
for ($n = 0;$n<=$len-25;$n++){
	my $temp_str = substr $consensus,$n,25;
	my $count = ($temp_str =~ tr/ATCG/ATCG/);
	if ($count >=22){
		$tag ++;
		$location{$tag} = $n;
	}
}
#3######################################################################
my %pairs;
my @tags = sort {$a<=>$b} keys %location;
if (@tags >= 2 ) {
	my ($A,$B);
	while(@tags) {
		my $A = shift @tags;
		foreach $B (@tags) {
#			print "$A\t$B\n";
			my $expected = $location{$B}-$location{$A};
			if ( $expected >= 425 && $expected <=825) { ##Here expected is PCR product and fwd primer.
				$pairs{$A}{$B} = $expected;
		#		print "$location{$B}\t$location{$A}\n";
			}
		}
	}
}
#4##############################
my ($key1,$key2);
my @sorted_array1 = sort{$a<=>$b}keys %pairs;
foreach $key1 (@sorted_array1) {
	my $hash2 = $pairs{$key1};
	my @sorted_array2 = sort{$a<=>$b}keys %$hash2;
		foreach $key2 (@sorted_array2) {
	#	print "$key1\t$key2\n";
			my @distances = &distance($location{$key1},$location{$key2});
			my $max = max(@distances);
			my $min = min(@distances);
			$max = sprintf "%.3f",$max;
			$min = sprintf "%.3f",$min;
			my $f = substr($consensus,$location{$key1},25);
			my $r = substr($consensus,$location{$key2},25);
			my $left_e = $location{$key1}+24;
			my $right_e = $location{$key2}+24;

			printf "#$location{$key1}-$left_e\t$f\t|\t$location{$key2}-$right_e\t$r\t$max\t$min\n";
		}
}
######################################################
sub hybird {
	my @a = @_;
	my $jianb ;
	if (($a[0] eq "A" && $a[1] eq "G") || ($a[1] eq "A" && $a[0] eq "G")) {
		$jianb = "R";
	}elsif (($a[0] eq "A" && $a[1] eq "C") || ($a[1] eq "A" && $a[0] eq "C")){
		$jianb = "M";
	}elsif (($a[0] eq "A" && $a[1] eq "T") || ($a[1] eq "A" && $a[0] eq "T")){
		$jianb = "W";
	}elsif (($a[0] eq "G" && $a[1] eq "C") || ($a[1] eq "G" && $a[0] eq "C")){
		$jianb = "S";
	}elsif (($a[0] eq "G" && $a[1] eq "T") || ($a[1] eq "G" && $a[0] eq "T")){
		$jianb = "K";
	}elsif (($a[0] eq "C" && $a[1] eq "T") || ($a[1] eq "C" && $a[0] eq "T")){
		$jianb = "Y";
	}
	return $jianb;
}
########################################################
sub distance {
	my $start= shift;
	my $end = shift;
	my @parts;
	foreach my $str(values %hash) {
		my $part = substr($str,$start+25,$end-$start-25);
		push @parts,$part;
	}
	my $que = shift @parts;
	my @distances;
	foreach my $obj (@parts) {
		my $diverge = &dis($que,$obj);
		push @distances,$diverge;
	}
	return @distances;
}
###############################################

sub dis {
	my @seq = @_;
	my $length = length ($seq[0]);
	my ($y,@array);
	my $conserved = 0;
	my $gap =0;
	for($y=0;$y<$length;$y++){
		foreach my $align(@seq){
			my $aa=substr $align,$y,1;
			push @array,$aa;
		}
		my %count;
		my @tmp=grep { ++$count{ $_ } < 2; } @array;
		if (@tmp==1 && $tmp[0]=~/[ATCG]/) {
			$conserved++;
		}elsif (grep /-/,@tmp){
			$gap++;
		}
		@array=();
	}
	my $diverg = ($length-$gap-$conserved)/($length-$gap);
	return $diverg;
}
