use List::Util qw(max);
###############sub##############
sub repeat_region{
my $lens=shift;
my $long=length ($lens);
my  $j=0;my @count;my $i;
for ($i=0;$i<=$long;$i++){
	my $base=substr $lens,$i,1;
	if ($base eq "*"){
	$j++;
	}else{
	push @count,$j;
	$j=0;
	}
}
my $max_repeat=max(@count);
return $max_repeat;
}
###################################
open IN,shift;
my $len=<IN>;
my $n=repeat_region($len);
print "$n";

