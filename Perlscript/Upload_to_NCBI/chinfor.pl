#!/usr/bin/perl 
use strict;
die "Usage: < part all>!\n"
unless @ARGV;
my @infor= &openfile($ARGV[0]);
#############################################################################
foreach(@infor){
	my @a=(split /\t+/,$_);
	if($a[0]=~/\#/){
		print;
		next;	
	}else{

		my $readlen=(split /\//,$a[6])[0];
		$a[11]=$a[11]/2/$readlen;	
#		print "$a[11]\t$readlen\n";
		print join "\t",@a; 

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
