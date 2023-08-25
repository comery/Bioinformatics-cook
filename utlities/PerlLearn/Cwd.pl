use Cwd qw(abs_path getcwd);
my $get = getcwd();
my $abs = abs_path("formatOutput.pl");
print "the abs path is : $abs\n";
print "the current path is : $get\n";
