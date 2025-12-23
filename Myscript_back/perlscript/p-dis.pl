#!usr/bin/perl -w 
###################################################
#this perl is to change the format of muscle output
#both DNA and PRO seq can be accepted!
###################################################
use strict;
use List::Util qw(min max);
open IN,$ARGV[0];
open OUT, ">>marker_distance.txt";
my $marker = (split /_/,$ARGV[0])[0].(split /_/,$ARGV[0])[1].(split /_/,$ARGV[0])[2];
my ($xingxing,%hash,$id,$seq);
$/=">";<IN>;$/="\n";
while ($id = <IN>){
	chomp $id;
	$/=">";
	$seq = <IN>;
	chomp $seq;
	$hash{$id} = $seq;
	$/="\n";
}

################### print blank format ############
my @head = keys %hash;
my (%distance,@distances);
while (@head) {
	my $que = shift (@head);
	foreach my $key (@head) {
	#	$distance{$que}{$key} = &dis($hash{$que},$hash{$key});
		my $diverg = &dis($hash{$que},$hash{$key});
		push @distances,$diverg;
#		print "$que <----> $key\t$diverg\n";
	}
}
my $min_dis = min (@distances);
my $max_dis = max (@distances);
print OUT "$marker\t$max_dis\t$min_dis\n";
close IN;
close OUT;
#####################################################
sub dis {
	my @seq = @_;
	my $length = length ($_[0]);
	my ($y,@array);
	my $conserved = 0;
	for($y=0;$y<$length;$y++){
		foreach my $align(@seq){
			my $aa=substr $align,$y,1;
			push @array,$aa;
		}
		my %count;
		my @tmp=grep { ++$count{ $_ } < 2; } @array;
		if (@tmp==1 && $tmp[0]=~/[A-Z]/) {
	#		$symble.="*";
			$conserved++;
		}
		@array=();
	}
	my $diverg = ($length-$conserved)/$length;
	return $diverg;
}
