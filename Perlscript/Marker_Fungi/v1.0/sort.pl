#!usr/bin/perl -w
use strict;
open IN,shift;
die "Usage:perl sort.pl id.list" unless (!@ARGV==1);
my (@sort,$id);
while ($id=<IN>){
chomp $id;
my $fam_id=$2 if ($id=~m/(\S+)_(\w+)$/);
push @sort,$id;
	if (@sort==5){
	print "$fam_id\t@sort\n";
	@sort=();
	}
}

