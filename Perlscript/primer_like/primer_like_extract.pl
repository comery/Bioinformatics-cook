#!/usr/bin/perl 

=head1 Name
	
	primer_like.pl - use a block to hit conserved regions of the sequences just
    like use a series of primers to extract expected DNA sequences.

=head1 Description
	
	Approximate string matching using dynamic programming
	Find the closest match to the pattern in the seq
	Gap will be considered more than mismatch.
	--p(primer) primer must shorter than sequences in --fa

=head1 Version
	
	Author: ChentaoYang, yangchentao@genomics.cn
	Version: 1.0,  Date:2015-1-28

=head1 Usage
	
	perl $0  [option]
	--p <str> primer string
	--fa <str>	subject string
	--g <num> gap punish
	--m <num> mismatch punish
	--cm <num> mismatch allowable(<=)
	--cg <num> gap allowable(<=)
	--o <str> output file including best hits
	--log <str>	log information including matched pattern and locations

=head1 Exmple

	perl ./primer_like_extract.pl --p primer.txt --fa seq.txt -o output_file -log log_file &
	perl ./primer_like_extract.pl --p primer.txt --fa seq.txt --gap_punish 4 --m 2 -cm 4 -cg 2 -o output_file -log log_file
	primer.txt may be like this: EIQADEVRL 
	seq.txt may be like this: 
	>id1
	SVLQDRSMPHQEILAADEVLQESEMRQQDMISHDE
	>id2
	SLQDRSMPHQEILABCDEVLQESEMRQQDMISHDE
=cut

use strict;
use Getopt::Long;

my ($gap_punishment,$mis_match,$primer,$fa,$cutoff_mis,$cutoff_gap,$log,$output,$help);
GetOptions(
	"p:s" => \$primer,
	"fa:s" => \$fa,
	"g:i" => \$gap_punishment,
	"m:i" => \$mis_match,
	"cm:i" => \$cutoff_mis,
	"cg:i" => \$cutoff_gap,
	"log:s" => \$log,
	"o:s" => \$output,
	"help!" => \$help
);

$cutoff_mis||= 2;
$cutoff_gap||= 1;
$gap_punishment ||= 2;
$mis_match||= 1;
$output||="$fa.out";
$log||="check.$fa.log";
die `pod2text $0` if (!($primer && $fa) || $help);
print "$cutoff_mis\t$cutoff_gap\n";
open QUE,"$primer";
open OBJ,"$fa";
open OUT,">$output"; 
open LOG,">$log";
my $pattern = <QUE>;
chomp $pattern;
my $PLEN = length $pattern;

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
	# D is the Distance matrix, which shows the "edit distance" between
	# substrings of the pattern and the seq.
	# It is implemented as a reference to an anonymous array.
	$D = [  ];

	# The rows correspond to the seq
	# Initialize row 0 of D.
	for (my $t=0; $t <= $TLEN ; ++$t) {
		$D->[$t][0] = 0;
	}

	# The columns correspond to the pattern
	# Initialize column 0 of D.
	for (my $p=0; $p <= $PLEN ; ++$p) {
		$D->[0][$p] = $p*$gap_punishment;
	#	$D->[0][$p] = 0;
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

	# Print D, the resulting edit distance array
	print LOG "//$id\n";	
#	for (my $p=0; $p<=$PLEN; $p++){
#		print LOG "        "if($p==0);
#		print LOG substr($pattern,$p-1,1), "   " if($p>0);
#	}
#	print LOG "\n";
#	for (my $t=0; $t <= $TLEN ; ++$t) {
#		print LOG substr($seq,$t-1,1), " " if($t>0);
#		print LOG "  "if($t==0);
#		for (my $p=0; $p <= $PLEN ; ++$p) {
#			printf LOG "%3.0f",$D->[$t][$p]; print LOG " ";
#		}
#		print LOG "\n";
#	}

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

	#print "PATTERN:\n$pattern\n";
	#print "TEXT:\n$seq\n";
	#print "gap_punish:$gap_punishment\n";
	#print "mismatch_punish:$mis_match\n";
	# Report the best match(es).
	#print "\nThe best match for the pattern $pattern\n";
	#print "has an edit distance of $bestscore\n";
	#print "and appears in the seq ending at location";
	#print "s" if ( @matches > 1);
	#print " @matches\n";


	# Print out
	my ($p_out, $t_out,%match_out);
	my $mark = 0;
	foreach my $matches (@matches){
		my $plen=$PLEN;my $matches_ori=$matches;
		$p_out='';$t_out='';
		$t_out=(substr ($seq, $matches-1, 1)).$t_out;
		$p_out=(substr ($pattern, $plen-1, 1)).$p_out;
		while($matches>1 && $plen>1){
			($seq,$pattern,$matches,$plen,$t_out,$p_out,$D)= &add_world ($seq,$pattern,$matches,$plen,$t_out,$p_out,$D);
		}
		print LOG "pattern_location:$plen-$PLEN	$p_out\n";
		print LOG "seq_location:$matches-$matches_ori	$t_out\n";	
		
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
		if (($mismatch_count < $cutoff_mis) && ($gap_count < $cutoff_gap)) {
			$mark++;
		}
	}
	if ($mark > 0){ print OUT ">$id\n$seq\n" ;$mark = 0 ;} #output the result which are up to standard.
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

