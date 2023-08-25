#!usr/bin/perl -w
my $Usage ="
	perl $0 <*_1.fq.gz> <*_2.fq.gz> <outprefix> <length|16bp>
";
use strict;

my $trim;
if (@ARGV < 4) {
	print "$Usage";
	exit();
}elsif (@ARGV == 3){
	$trim = 16;
}else{
	$trim = $ARGV[3];
}
my $prefix = $ARGV[2];
#input
if ($ARGV[0] =~ /.gz$/) {
	open IN1,"gzip -dc $ARGV[0]|" or die "Can't open $ARGV[0] !";
}else {
	open IN1,"$ARGV[0]" or die "Can't open $ARGV[0] !";
}
if ($ARGV[1] =~ /.gz$/) {
	open IN2,"gzip -dc $ARGV[1]|" or die "Can't open $ARGV[1] !";
}else {
	open IN2,"$ARGV[1]" or die "Can't open $ARGV[1] !";
}

#output
open OUT1,"|gzip > $ARGV[2]\_1.fq.gz" or die "Can't open $ARGV[2]\_1.fq.gz !";
open OUT2,"|gzip > $ARGV[2]\_2.fq.gz" or die "Can't open $ARGV[2]\_2.fq.gz !";

my $count = 0;
my ($line1,$line2);

while ($line1 = <IN1> ) {
	$line2 = <IN2>;
	$count ++;
	my $id1 = &order($line1);
	my $id2 = &order($line2);
	if ($id1 ne $id2) {
		die "fq1 and fq2 are at worry order! Check it out!";
	}
	my $seq1 = <IN1>;
	my $seq2 = <IN2>;
	# +
	<IN1>;<IN2>;
	my $qual1 = <IN1>;
	my $qual2 = <IN2>;
	my $adapter = substr($seq1,0,$trim);
	substr($seq1,0,$trim) = "";
	substr($qual1,0,$trim) = "";
	
	print OUT1 "\@$prefix-$count\_$adapter\n";
	print OUT1 "$seq1";
	print OUT1 "+\n$qual1";
	print OUT2 "\@$prefix-$count\_$adapter\n";
	print OUT2 "$seq2";
	print OUT2 "+\n$qual2";
}

close IN1;
close IN2;
close OUT1;
close OUT2;


sub order {
	my $str1 = shift;
	my $lane1 = (split /\s+/,$str1)[0];
	return $lane1;	
}
