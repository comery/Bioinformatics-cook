#!/usr/bin/perl 

=head1 Name
	
	./primer_like_extract.pl for HIFI barcoding

	 - use a block to hit conserved regions of the sequences just
    like use a series of primers to extract expected DNA sequences.

=head1 Description
	subdivide CCS(circular consensus sequence) to sample by index which is 
	linked primers,like below:
    type 1:
    [AAAGC]GGTCAACAAATCATAAAGATATTGG------------TGATTTTTTGGTCACCCTGAAGTTTA[GCTTT]

    type 2:
    [AAAGC]TAAACTTCAGGGTGACCAAAAAATCA------------CCAATATCTTTATGATTTGTTGACC[GCTTT]

	 here [AAAGC] is index for one sample, while [GCTTT] is reverse and complementary.

	Approximate string matching using dynamic programming
	Find the closest primer matches in the seq
	Gap will be considered more than mismatch.
	In addition, index regions don't allow mismatch or gap,
	just considering 1 bp loss in head or tail.
	--p(primer) primer must shorter than sequences in --fa

=head1 Version
	
	Author: ChentaoYang, yangchentao@genomics.cn
	Version: 1.0,  Date:2015-1-28
	Version: 2.0,  Date:2017-06-01

=head1 Usage
	
	perl $0  [option]
	--p <str> primer string
	--index <str> samples' index sequence
	--fa <str>	subject string
	--g <num> gap punish
	--m <num> mismatch punish
	--cm <num> mismatch allowable(<=)
	--cg <num> gap allowable(<=)
	--o <str> output file including best hits
	--outdir <outdir> outdir contains subdivided samples' ccs
	--log <str>	log information including matched pattern and locations

=head1 Exmple

	perl ./primer_like_extract.pl --p primer.txt --index index.xls --fa ccs.fa -o output_file  -log log_file &
	perl ./primer_like_extract.pl --p primer.txt --index index.xls --fa ccs.fa --gap_punish 4 --m 2 -cm 4 -cg 2 -o output_file  -log log_file
	
	primer.txt may be like this:
	
		for	GGTCAACAAATCATAAAGATATTGG
		rev	TAAACTTCAGGGTGACCAAAAAATCA

	!!! ATTENTION:!!!!!!!!!!!!!!!!!!!!!!
		for-revcom(rev) is a pari regarding to type 1
		rev-revcom(for) is a pari regarding to type 2
	
	ccs.fa may be like this: 
	
	>m170506_092957_42199_c101149142550000001823255607191735_s1_p0/25/ccs
	GTTACTAAACTTCAGGGTGACCAAAAATCAAAATAGATGTTGATAGAGGATAGGGTCACCTCCTCCGCAG
	GGTCAAAAATGAAGTATTCAAATTTCGATCTGTTAAAAGTATTGTAATAGCTCCGGCTAAAACTGGTAAT
	GAGAGAAGGAGTAAAAGGCAGTAATGACGACTGATCAGACAAATAGGGTATTCGATCAAGGGTAATCCGG
	GAGATCGTATGTTAATCACTGTTGTGATGAAATTTACGGCCCTAAAATAGAGGAAATTCCGGCTAAATGA
	AGGGAAAAAATAGCTAAGTCAACAGAAGCCCCTCTATGGGCAATTCCCGAAGAAAGGGGGGATAAACCGT
	TCACCCAGTTCCTGCCCCATTTCAACTATACTTCTAGCTAAAAGAAGGATAAAGAAGGGGGAGTATTCAA
	AATCTTATATTATTCATTCGAGGAAAGGCCATGTCTGGTGCTCCTAATATGAGAGGGACAAGTCAGTTTC
	CAAATCCCCAATTATAATAGGTATAACCATAAAAAAAATTATAATAAAAGCGTGGCTGTAACAATTACAT
	TATAGATTTGATCATCTCCAATGAGAGAGCCGGGATGTCCTAATTCAGTACGAATTAATATACTTAAAGA
	AGTACCAACTATCCAGCCCATGCCCCAAAGAGAAAGTAAAGAGTACCAATATCTTTATGATTTGTTGACC
	GTAA

=cut

use strict;
use Getopt::Long;

my ($gap_punishment,$mis_match,$indexlst,$primer,$fa,$cutoff_mis,$cutoff_gap,$log,$output,$help);
GetOptions(
	"p=s" => \$primer,
	"index=s" => \$indexlst,
	"fa=s" => \$fa,
	"g:i" => \$gap_punishment,
	"m:i" => \$mis_match,
	"cm:i" => \$cutoff_mis,
	"cg:i" => \$cutoff_gap,
	"log:s" => \$log,
	"o:s" => \$output,
	"help!" => \$help
);
$cutoff_mis ||= 2;
$cutoff_gap ||= 1;
$gap_punishment ||= 2;
$mis_match ||= 1;
$output ||="ccs.success.fa";
$log ||="check.$fa.log";
die `pod2text $0` if (!($primer && $fa && $indexlst ) || $help);
die "So embarrassing! cutoff_mis will be set to 2, and cutoff_mis will be 1 in spite you want 0" if ($cutoff_mis == 0 && $cutoff_gap ==0);

open QUE,"$primer";
open OBJ,"$fa";
open OUT,">$output"; 
open LOG,">$log";
open LOG1,">split.log.txt";

print LOG "###id\tprimer_match_type\tmismatch_count\tgap_count\n";
my %patt;
while (<QUE>) {
	chomp;
	my @a = split;
	$patt{'1'} = $a[1] if ($a[0] eq "for") ;
	$patt{'2'} = $a[1] if ($a[0] eq "rev") ;
}
my $patt_n = keys %patt;
unless ($patt_n == 2){
	print  "There are more than two lines in primer list!\nHere is example:\n
	\tfor GGTCAACAAATCATAAAGATATTGG
	\trev TAAACTTCAGGGTGACCAAAAAATCA";
	exit;
}


###add 3,4 into %patt;
$patt{'3'} = comrev($patt{'2'});
$patt{'4'} = comrev($patt{'1'});
#----------------------------------------------------------------
open INDEX,$indexlst;
my (%indexf,%indexr,%trimedf,%trimedr);
while (<INDEX>) {
	chomp;
	next if (/^\s+$/);
	my @a = split /\s+/,$_;
	my $num = $a[0];
	my $pri = $a[1];
	$indexf{$pri} = $num;
#	print "5 nor\t$pri\t$num\n";
	my $trim_pri = &trim_head($pri);
	$trimedf{$trim_pri} = $num;
#	print "4 nor\t$trim_pri\t$num\n";
	my $revcom = &comrev($pri);
	$indexr{$revcom} = $num;
#	print "5 rev\t$revcom\t$num\n";
	my $trim_revcom = &trim_tail($revcom);
	$trimedr{$trim_revcom} = $num;
#	print "4 rev\t$trim_revcom\t$num\n";
	
}
close INDEX;

#Loop each sequence for check
my $D;
$/=">";<OBJ>;$/="\n";
while (my $id = <OBJ>) {
	chomp $id;
	$/=">";
	my $seq = <OBJ>;
	chomp $seq;
	$seq =~ s/\s//g;
	$seq =~ s/\n//g;
	my $TLEN = length $seq;
	my (%mark , %border_arr,%border);
	foreach my $key(keys %patt){
		my $pattern = $patt{$key};
		my $PLEN = length $pattern;
		$D = [  ];

		for (my $t=0; $t <= $TLEN ; ++$t) {
			$D->[$t][0] = 0;
		}

		for (my $p=0; $p <= $PLEN ; ++$p) {
			$D->[0][$p] = $p*$gap_punishment;
		}

		# Compute the edit distances.
		for (my $t=1; $t <= $TLEN ; ++$t) {
			for (my $p=1; $p <= $PLEN ; ++$p) {

				$D->[$t][$p] =

		# Choose whichever of the three alternatives has the least cost
				min3(
		# First alternative
		# The seq and pattern may or may not match at this character ...
						substr($seq, $t-1, 1) eq substr($pattern, $p-1, 1) 
						? $D->[$t-1][$p-1]  # If they match, no increase in edit distance!
						:  $D->[$t-1][$p-1] + $mis_match,

		# gap +1 mismatch +1 
		# Second alternative
		# If the seq is missing a character
						$D->[$t-1][$p] + $gap_punishment, 

		# Third alternative
		# If the pattern is missing a character
						$D->[$t][$p-1] + $gap_punishment
				    )
			}
		}


		my @matches = (  );
		my $bestscore = 10000000;

		# Find the best match(es).
		# The edit distances appear in the the last row.)
		for (my $t=1 ; $t <= $TLEN ; ++$t) {
			if( $D->[$t][$PLEN] < $bestscore) {
				$bestscore = $D->[$t][$PLEN];
				@matches = ($t);
			}elsif( $D->[$t][$PLEN] == $bestscore) {
				push(@matches, $t);
			}
		}

		# Print out
		my ($p_out, $t_out,%match_out);
		
		foreach my $matches (@matches){
			my $plen=$PLEN;my $matches_ori=$matches;
			$p_out='';$t_out='';
			$t_out=(substr ($seq, $matches-1, 1)).$t_out;
			$p_out=(substr ($pattern, $plen-1, 1)).$p_out;
			while($matches>1 && $plen>1){
				($seq,$pattern,$matches,$plen,$t_out,$p_out,$D)= &add_world ($seq,$pattern,$matches,$plen,$t_out,$p_out,$D);
			}
			
			
			# Check whether this sequence is ok
			my $len_out = length($t_out);
			my $mismatch_count = 0;my $gap_count = 0;
			for (my $x = 0,my $y = 0;$x<=$len_out-1,$y<=$len_out-1;$x++,$y++) {
				my $tmp_x = substr($p_out,$x,1);
				my $tmp_y = substr($t_out,$y,1);
				$mismatch_count++ 	if (($tmp_x ne '_') && ($tmp_y ne '_') && ($tmp_x ne $tmp_y));
				$gap_count++ 		if ($tmp_y eq '_' || $tmp_x eq '_');
			}
			#print $mismatch_count,"\n",$gap_count;

			# how many matched regions you can get when you set cutoff of mismatch and gap
			# $len_out == $PLEN is modified by yangchentao,20170627, 
			# because sometimes there is short match than patterm length, and that is wrong!
			my $real_len_out = $PLEN - $plen + 1;
			print "$real_len_out\t$PLEN\t\n";
			if (($mismatch_count <= $cutoff_mis) && ($gap_count <= $cutoff_gap) && ($real_len_out == $PLEN)) {
				$mark{$key}++;
				print LOG "//$id\t$key\t$mismatch_count\t$gap_count\n";
				print LOG "pattern_location:$plen-$PLEN	$p_out\n";
				print LOG "seq_location:$matches-$matches_ori	$t_out\n";	
				push @{$border_arr{$key}},$matches if ($key == 1 or $key == 2); #head part
				push @{$border_arr{$key}},$matches_ori if ($key == 3 or $key == 4); # tail part
			}else{
				print  "//$id\t$key\t$mismatch_count\t$gap_count\n";
				print  "pattern_location:$plen-$PLEN	$p_out\n";
				print  "seq_location:$matches-$matches_ori	$t_out\n";
			}
		}
		foreach my $j(sort keys %border_arr){
		
			if ($j == 1 or $j == 2){my $min = $border_arr{$j}->[0] ;$border{$j} = $min;}

			if ($j == 3 or $j == 4){my $max = $border_arr{$j}->[-1] ;$border{$j} = $max;}
			
		}
		
	}

#	print "$id\t$mark{'1'}\t$mark{'3'},$mark{'2'}\t$mark{'4'}\n";
	my ($remanet,$head,$tail);
	if ($mark{'1'} * $mark{'3'} >0 ) { # $mark{'1'} >0 && $mark{'3'} >0, means that there are at least one for-primer match and one rev-primer match
		my $head_len = $border{'1'}-1;
		$remanet = $TLEN - $border{'3'};
		$head = substr($seq,0,$head_len);
		$tail = substr($seq,$border{'3'},$remanet);
		if ($head_len < 4 || $remanet < 4 ) {
			print LOG1 "$id can't be belonged to any sample! because of short index[<4bp]\n";
		}else {
			my $num1 = &own1($head,$head_len);
			my $num2 = &own2($tail,$remanet);
			print LOG1 "$id\t$head\t$head_len\t$num1\t$tail\t$remanet\t$num2\n";
		#	open OUT,"$outdir/$num1/$num1.for.fa";
			print OUT ">$id\t$num1\tfor\n$seq\n" if ($num1 == $num2 && $num1 != 0);
		#	close OUT,
		}
		
	#	print OUT ">$id\n$seq\n" ;
	}elsif($mark{'2'} * $mark{'4'} > 0) { # $mark{'2'} >0 && $mark{'4'} >0
		my $head_len = $border{'2'}-1;
		$remanet = $TLEN - $border{'4'};
		$head = substr($seq,0,$head_len);
		$tail = substr($seq,$border{'4'},$remanet);
		if ($head_len< 4 || $remanet < 4 ) {
			print LOG1 "$id can't be belonged to any sample! because of short index[<4bp]\n";
		}else {

			my $num1 = &own3($head,$head_len);
			my $num2 = &own4($tail,$remanet);
			print LOG1 "$id\t$head\t$head_len\t$num1\t$tail\t$remanet\t$num2\n";
			print OUT ">$id\t$num1\trev\n$seq\n" if ($num1 == $num2 && $num1 != 0 );
		}
	}elsif($mark{'1'} + $mark{'3'} > 2 || $mark{'2'} + $mark{'4'} > 2) {
		print LOG1 "$id has too many primer regions!\n"; ## error, more than one for primer or rev primer
	}else { ## no matches sastifying your cutoff
		print LOG1 "$id no primer region under your cutoff\n";
	}

	%border = ();
	%mark = ();
	$/="\n";

}

	

sub min3 {
	my($i, $j, $k) = @_;
	my($tmp);

	$tmp = ($i < $j ? $i : $j);
	$tmp < $k ? $tmp : $k;
}
	
sub add_world {
	my ($seq_s,$pattern_s,$matches_s,$PLEN_s,$t_out_s,$p_out_s,$D_s)=@_;
	my $flag=0;
# $flag=1 <-      2 |     3 \
	if(substr ($pattern_s, $PLEN_s, 1) ne substr ($seq_s, $matches_s, 1)){
		if($flag==0 && $D_s->[$matches_s][$PLEN_s-1]+$gap_punishment <= $D_s->[$matches_s-1][$PLEN_s-1]+$mis_match && $D_s->[$matches_s][$PLEN_s-1] <= $D_s->[$matches_s-1][$PLEN_s]){
			$PLEN_s=$PLEN_s-1;
			$t_out_s="_".$t_out_s;
			$p_out_s=(substr ($pattern_s, $PLEN_s-1, 1)).$p_out_s;
			$flag=1;
		}


		if($flag==0 && $D_s->[$matches_s-1][$PLEN_s] <= $D_s->[$matches_s][$PLEN_s-1] && $D->[$matches_s-1][$PLEN_s]+$gap_punishment <= $D_s->[$matches_s-1][$PLEN_s-1]+$mis_match){
			$matches_s=$matches_s-1;
			$t_out_s=(substr ($seq_s, $matches_s-1, 1)).$t_out_s;
			$p_out_s="_".$p_out_s;
			$flag=2;
		}
		if($flag==0 && $D_s->[$matches_s-1][$PLEN_s-1]+$mis_match <= $D_s->[$matches_s][$PLEN_s-1]+$gap_punishment && $D_s->[$matches_s-1][$PLEN_s-1]+$mis_match <= $D_s->[$matches_s-1][$PLEN_s]+$gap_punishment){
			$matches_s=$matches_s-1;
			$PLEN_s=$PLEN_s-1;
			$t_out_s=(substr ($seq_s, $matches_s-1, 1)).$t_out_s;
			$p_out_s=(substr ($pattern_s, $PLEN_s-1, 1)).$p_out_s;
			$flag=3;
		}
	}
	if(substr ($pattern_s, $PLEN_s, 1) eq substr ($seq_s, $matches_s, 1)){
		if($flag==0 && $D_s->[$matches_s][$PLEN_s-1]+$gap_punishment <= $D_s->[$matches_s-1][$PLEN_s-1] && $D_s->[$matches_s][$PLEN_s-1] <= $D_s->[$matches_s-1][$PLEN_s]){
			$PLEN_s=$PLEN_s-1;
			$t_out_s="_".$t_out_s;
			$p_out_s=(substr ($pattern_s, $PLEN_s-1, 1)).$p_out_s;
			$flag=1;
		}
		if($flag==0 && $D_s->[$matches_s-1][$PLEN_s] <= $D_s->[$matches_s][$PLEN_s-1] && $D->[$matches_s-1][$PLEN_s]+$gap_punishment <= $D_s->[$matches_s-1][$PLEN_s-1]){
			$matches_s=$matches_s-1;
			$t_out_s=(substr ($seq_s, $matches_s-1, 1)).$t_out_s;
			$p_out_s="_".$p_out_s;
			$flag=2;
		}
		if($flag==0 && $D_s->[$matches_s-1][$PLEN_s-1] <= $D_s->[$matches_s][$PLEN_s-1]+$gap_punishment && $D_s->[$matches_s-1][$PLEN_s-1] <= $D_s->[$matches_s-1][$PLEN_s]+$gap_punishment){
			$matches_s=$matches_s-1;
			$PLEN_s=$PLEN_s-1;
			$t_out_s=(substr ($seq_s, $matches_s-1, 1)).$t_out_s;
			$p_out_s=(substr ($pattern_s, $PLEN_s-1, 1)).$p_out_s;
			$flag=3;
		}
	}
	return ($seq_s,$pattern_s,$matches_s,$PLEN_s,$t_out_s,$p_out_s,$D_s) if($flag!=0);
}


sub comrev {
	my $a = shift;
#	print "$a\n";
	chomp $a;
	$a =~ tr/NATCG/NTAGC/;
	$a= reverse $a;
	return $a;
}

sub trim_tail {
	my $b = shift;
	$b =~ s/(\w)$//;
	return $b;
}

sub trim_head {
	my $c = shift;
	$c =~ s/^(\w)//;
	return $c;
}

sub own1 {
	my ($str,$len) = @_;
	if ($len >= 5) {
		my $tmp = substr($str,-5);
		exists $indexf{$tmp} ? return $indexf{$tmp} : return 0;
	}else {
		exists $trimedf{$str} ? return $trimedf{$str} : return 0;
	}
}

sub own2 {
	my ($str,$len) = @_;
	if ($len >= 5) {
		my $tmp = substr($str,0,5);
		exists $indexr{$tmp} ? return $indexr{$tmp} : return 0;
	}else {
		exists $trimedr{$str} ? return $trimedr{$str} : return 0;
	}
}

sub own3 {
	my ($str,$len) = @_;
	if ($len >= 5) {
		my $tmp = substr($str,-5);
		exists $indexf{$tmp} ? return $indexf{$tmp} : return 0;
	}else {
		exists $trimedf{$str} ? return $trimedf{$str} : return 0;
	}
}

sub own4 {
	my ($str,$len) = @_;
	if ($len >= 5) {
		my $tmp = substr($str,0,5);
		exists $indexr{$tmp} ? return $indexr{$tmp} : return 0;
	}else {
		exists $trimedr{$str} ? return $trimedr{$str} : return 0;
	}
}

