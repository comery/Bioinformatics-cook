#!/usr/bin/perl 
use strict;
die "Usage: < part all>!\n"
unless @ARGV >=2;
my @infor= &openfile($ARGV[0]);
my %hash=&index;
my %stat;
#############################################################################
foreach(@infor){
	chomp;
	my ($id,$gene)=(split /\s+/,$_)[1,0];
	$id=($gene=~/\d+_([^_]+)_\d+/)[0];

	if(exists $hash{$id}){
		print "$hash{$id}\t$_\n"; 
		$stat{$id}++;
	}
}
foreach my $id(keys %stat){
	if($stat{$id} == 1){
		print "NA\t$hash{$id}\n";
	}


}
##############################################################################
sub index{
	my @go = &openfile($ARGV[1]);
	my %hash;
	foreach(@go){
		chomp;
		my($id,$infor)=(split /\s+/,)[1,1];
		$hash{$id}=$_; # if $p > 0.5;
		$stat{$id}=1;
	}
	return %hash;
}
sub openfile{
	my ($filename)=@_;
	open LISTFILE ,$filename or die "you can not open $filename\n";
	my @seq=<LISTFILE>;
	close LISTFILE ;
	return @seq;
}
