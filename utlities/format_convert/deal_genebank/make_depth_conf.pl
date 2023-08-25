use strict;
die "Usage : perl $0 <*.depth> >*.depth.txt" unless (@ARGV ==1);
open IN,shift;
while (<IN>) {
	next if (/average_depth/);
	chomp;
	my @a = split /\s+/,$_;

	foreach my $i (3..$#a){

		my $j = $i - 2 ;
		my $n = $j + 1;
		print "mt1\t$j\t$j\t$a[$n]\n";
	}
}

