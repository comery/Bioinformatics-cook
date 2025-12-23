#!/usr/bin/perl -w
use strict;
if (@ARGV <2) {
	print "Usage:\n\tperl $0 <md5 file on linux> <md5 file on Mac>";
	exit();
}

open LIN,shift;
open MAC,shift;
my @arr;
my (%linux,%mac);

while (<LIN>) {
	chomp;
	my @a = split;
	$linux{$a[1]} = $a[0];
	push @arr,$a[1];
}
my $total = @arr;
while (<MAC>) {
	chomp;
	my @b = split /=/;
	my $file = $1 if ($b[0] =~ /\((.+)\)/);
	#print "$file\n";
	$mac{$file} = $b[1];
}

my $s = 1;
for my $i (@arr){
	if (exists $linux{$i} && exists $mac{$i}){
		print "$i\tOK" if ($linux{$i} eq $mac{$i});
	
	}else{
		print "files in two md5 file are not identical!\n";
		$s = 1;
	}
}
if ($s) {
	print "All $total files are ok!";
}
