#!usr/bin/perl -w 
###################################################
#this perl is to change the format of muscle output
#both DNA and PRO seq can be accepted!
###################################################
use strict;
use List::Util qw(min max);
my $file = shift;
open IN,"$file" or die "$!!";
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

#1####################################################
my @seq = values %hash;
my $length=length($seq[0]);
my ($y,@array);
my $consensus ;
for($y=0;$y<=$length-1;$y++){
	foreach my $align(@seq){
		my $aa=substr $align,$y,1;
		push @array,$aa;
	}
	my (%count,$spot);
	my @tmp=grep { ++$count{ $_ } < 2; } @array;
	if (@tmp==1 && $tmp[0]=~/[ATCG]/) {
		$spot = "$tmp[0]";
	}elsif (grep /N/,@tmp){
		$spot="-";
	}elsif (grep /-/,@tmp) {
		$spot="-";
	}elsif (@tmp==2 ){
		$spot= &hybird(@tmp);
	}elsif (@tmp==3 ){
		$spot = &hybird1(@tmp);
	}elsif (@tmp==4 ){
		$spot="N";
	}
	$consensus.=$spot;
	@array=();
}

#2###############################################################
my (%location,@standard);
my $n;
my $len = length ($consensus);
my $tag =0;
for ($n = 0;$n<=$len-21;$n++){
	my $temp_str = substr $consensus,$n,21;
	my $snp = substr $temp_str,10,1;
	my $snp_location = $n+10;
	my $count = ($temp_str =~ tr/ATCG/ATCG/);
	if ($count ==20 && $snp ne "A" && $snp ne "T" && $snp ne "C" && $snp ne "G" ){
		push @standard,$snp_location;
		print "$file\t$snp_location\t$temp_str\n";
	}
}

#3######################################################################
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
	}else {
		$jianb = "-";
	}
	return $jianb;
}
######################################################################
sub hybird1 {
	my @b = @_;
	my $degenerated ;
	if (!grep /A/,@b){
		$degenerated = "B";
	}elsif (!grep /G/,@b){
		$degenerated = "H";
	}elsif (!grep /C/,@b) {
		$degenerated = "D";
	}elsif (!grep /T/,@b) {
		$degenerated = "V";
	}else {
		$degenerated = "-";
	}
	return $degenerated;
}
#######################################################################

