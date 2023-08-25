#!/usr/bin/perl
=head1 Name
	deal_newick.pl

=head1 Contact & Version
		Author: Chentao YANG, yangchentao@genomics.org.cn
		Version: 1.0,  Date: 2015-11-24

=head1 Command-line Option
		--leaf <string>		leaf name of tree
		--clade <string>	a range clade(e.g.:A-B)
		--help		output help information to screen

=head1 Usage Exmples
		perl  ./deal_newick.pl  -leaf Acerentomon_sp  -clade Machilis_hrabei-Atelura_formicaria test.new
	
=head1 Descripiton
		This script is for handling the newick format file, especially deleting some leaves in the tree 
		without changing the topological strcture. You can feed this script with a parameter "-leaf" 
		following one species name or use "-clade" parameter following a name range just like 
		"Machilis_hrabei-Atelura_formicaria".

=cut
#################################################################################
use strict;
use warnings;
use Getopt::Long;
my ($leaf,$clade,$Help);
GetOptions(
	"leaf:s" => \$leaf,
	"clade:s" => \$clade,
	"help:s" => \$Help);
#die `pod2text $0` if (@ARGV==0);
#die `pod2text $0` if ($Help);

open IN ,shift or die "Can't open the file!";
my $file = <IN>;
my $str = $file;
$str =~ s/\(//g;
$str =~ s/\)//g;
my @names = split /\,/,$str;
my @deal;

#print $file;
#---------------------------check name in file----------------------------#
my ($clade_a,$clade_b);
if ($clade) {
	$clade_a = (split /-/,$clade)[0];
	$clade_b = (split /-/,$clade)[1];
#	print "$clade_a\t$clade_b\n";
	die  "Can't find this name in your file!" unless ($file =~ /$clade_a/ && $file =~ /$clade_b/);
	my $suffix_b = &suffix($clade_b);
	my $suffix_a = &suffix($clade_a);
	for (my $n=$suffix_a;$n<=$suffix_b;$n++) {push @deal,$names[$n];}
}
if ($leaf) {
	my @aa = split /,/,$leaf;
	foreach my $l (@aa) {
		die "Can't find this name in your file!" unless ($file =~ /$l/);
		push @deal,$l;
	}
}
print "@deal\n";
die "Which one you want to delete?" if (!$leaf && !$clade);
#-------------------------------------------------------------------------#

foreach my $sp (@deal) {
	
	my ($L,$R) =split /$sp/,$file ;
	my @left =split //,$L;
	my @right =split //,$R;
#	print "@left\n";
	$left[-1]="" if ($L =~ /,$/);
#	print "@left\n";
	my $NUM=0;
	for (my $x=($#left);$x>=0;$x--){
		if ($left[$x] =~ /\)/) {
			$NUM++;
		}
		if ($left[$x] =~ /\(/ && $NUM!=0){
			$NUM--;
		}
		if ($left[$x] =~ /\(/ && $NUM==0){
			$left[$x]=""; 
			last;
		}
	}
	$right[0]="" if ($R =~ /^,/);
	my $num=0;
	for (my $y=0;$y<=$#right;$y++){
		if ($right[$y] =~ /\(/){
			$num++;
		}
		if ($right[$y] =~ /\)/ && $num!=0){
			$num--;
		}
		if ($right[$y] =~ /\)/ && $num==0){
			$right[$y]="";
			last;
		}
	}
	$file=join("",@left[0..$#left],@right[0..$#right]);
}
print "$file" ; 



sub suffix {
	my $key = shift;
	my %hash ;
	for (my $j = 0;$j <= $#names;$j++) {
		$hash{$names[$j]} = $j;
	}
	my $suffix_num = $hash{$key};
	return $suffix_num;
}
