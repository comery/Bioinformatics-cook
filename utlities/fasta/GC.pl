#！/usr/bin/perl -w
use strict;
use Getopt::Long;

=head1 Version
	
	Author: ChentaoYang, yangchentao@genomics.cn
	Version: 1.0,  Date:2017-3-10

=head1 Usage
	
	perl $0  *.clw [option]
	--link	to link all contigs into one ?
	--l	whether to link all contigs into one by length(from large to small) ?
	--win	<number>	how long sequence calculate a GC contain (defult 100)
	--out	<str>	outputfile
	--help	show this help information
	
=head1 Exmple

	perl -link -l -win 1000 test.fa

=head1 Result
	


=cut

my ($link,$bylen,$win,$output,$help);
GetOptions(
	"link!"=>\$link,
	"l!"=>\$bylen,
	"win:s"=>\$win,
    "out:s"=>\$output,
    "help!"=>\$help
    );
die `pod2text $0` if ($help || @ARGV==0);
my $fa = $ARGV[0];
$output ||= "$ARGV[0].gc.txt";
$win ||= 100;

if (! $link ) {
	# body...
	if ($bylen) {
		# body...
		&mode2($fa);
	}else {
		&mode1($fa);
	}
}else {
	if ($bylen) {
		# body...
		&mode2($fa);
	}else{
		&mode3($fa);
	}
}
#mode1	no link and by sequence of its like
#mode2	link and by sequence of length
#mode3	link but not by sequence of length


sub mode1{
	my $file = shift;
	open IN, "$file" or die "Can not open $file!";
	open OUT ,">$output";
	$/=">";<IN>;$/="\n";
	while (my $id=<IN>) {
		chomp $id;
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		my $len = length($seq);
		for (my $i = 0; $i < $len - $win -1; $i+=$win) {
			my $tmp = substr($seq,$i,$win);
		#	print "$tmp\n";
			my $GC = &gc($tmp);
		#	print "$GC\n";
			my $start = $i +1 ;
			my $end = $i + $win ;
			print OUT "$id\t$start\t$end\t$GC\n";
		}
		$/="\n";
	}
	close IN;
	close OUT;
 }

sub mode2{

}

sub mode3{

}

sub gc{
	my $str = shift;
	my $count_g = $str =~ s/G/G/gi;
	my $count_c = $str =~ s/C/C/gi;
	my $str_len = length($str);
	my $contain = ($count_g+$count_c)/$str_len;
	return $contain;
}

