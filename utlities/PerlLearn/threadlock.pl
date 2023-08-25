#!/usr/bin/perl
use strict;
use threads;
use threads::shared;
use Data::Dumper;

my $val     : shared; # 共享变量
my %hash    : shared; # 共享数组
my @array   : shared; # 共享哈希

my $t1 = threads->create(\&test1);
my $t2 = threads->create(\&test2);

$t1->join; # 回收 t1 的线程
$t2->join;

print Dumper(\$val);
print Dumper(\@array);
print Dumper(\%hash);

sub test1 {
	lock ($val); lock (@array); lock (%hash);
	for ( 1 .. 1000 ){
		$val++;
		$array[0]++;
		$hash{test}++;
	}
}

sub test2 {
	lock ($val); lock (@array); lock (%hash);
	for ( 1 .. 1000 ){
		$val++;
		$array[0]++;
		$hash{test}++;
	}
}
