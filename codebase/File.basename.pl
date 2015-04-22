#!usr/bin/perl -w
use Cwd;
use File::Basename;
my $fullname="./perlscript";
my ($name,$path,$suffix,@suffixlist,$basename,$dirname);
$basename = basename($fullname); 
$dirname = dirname($fullname); 
print "$dirname\n$basename";
