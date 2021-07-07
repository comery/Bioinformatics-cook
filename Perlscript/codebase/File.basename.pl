#!usr/bin/perl -w
use Cwd;
use File::Basename;
my $fullname="/ifs4/NGB_ENV/USER/yangchentao/software/perlscript/codebase/File_Basename.pl";
my ($name,$path,$suffix,@suffixlist,$basename,$dirname);
$basename = basename($fullname); 
$dirname = dirname($fullname); 
$basename =~ s/\.\w+$//g;
print "$dirname\n$basename";
