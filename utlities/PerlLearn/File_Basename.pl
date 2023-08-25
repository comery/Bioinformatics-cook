#!/usr/bin/perl
use File::Basename qw(basename);

my($name, $path, $suffix)=File::Basename::fileparse(formatOutput.pl);
my $filename = basename("/ifs4/NGB_ENV/USER/yangchentao/software/perlscript/codebase/formatOutput.pl",  ".pl");

print "$name\n$path\n$suffix\n";
print "$filename";
