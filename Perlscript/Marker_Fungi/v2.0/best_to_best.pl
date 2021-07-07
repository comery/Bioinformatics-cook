#!usr/bin/perl -w 
use strict;
die "Usage: perl $0 <hmm.all.results> <blastp.best.out>" unless (@ARGV==2);
open IN,shift;
open IN1,shift;
my (@a,@b,%hash1,%hash2,$line);
while (<IN>) {
	chomp;
	@a=split /\s+/,$_;
	my $gene = (split/_/,$a[1])[0];
	$hash1{$a[0]}=$gene;
}
while (<IN1>){
	chomp;
	@b=split /\s+/,$_;
	$hash2{$b[0]}=$b[1];
}
my @k1=keys %hash1;
my @k2=keys %hash2;
my %a = map{$_ => 1} @k1;
my %b = map{$_ => 1} @k2;
my @inter = grep {$a{$_}} @k2;
#print "@inter\n";
foreach my $key(@inter){
		if($hash1{$key} eq $hash2{$key}){
		print "$key\t$hash1{$key}\n";
		} 
	}
close IN;
close IN1;
