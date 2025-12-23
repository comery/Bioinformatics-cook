#!usr/bin/perl -w 
use strict;
use List::Util qw(max);
#################################################################
# Get the consensu sequence for this segment
#################################################################

my $usage = " Usage:  perl $0 seg1-out.base.txt reference_cds_location.fa \n";

#######################################
#
# Test for commandline arguments
#
#######################################
if (scalar @ARGV != 2) {
 	print "$usage";
	exit -1;
}

my $file = $ARGV[0];
open IN,"$file";##<input file> Influenza-A-seg.out
my (@info,$GI,%sites,@sub_sites,@lines,$line,$len,$key,%key_base,$sequence);
while (<IN>){
	next if ($_ =~ /Scaffold/ );
	@info=split /\s+/,$_;
	$GI = $1 if ($info[0] =~ /gi\|(\d+)\|gb/);
#	print "$GI\n";
	$key=$info[1];
	my $aa="A";my $t="T";my $g="G";my $cc="C";
	push @sub_sites,$aa if ($info[3]>0);
	push @sub_sites,$t if ($info[4]>0);
	push @sub_sites,$g if ($info[5]>0);
	push @sub_sites,$cc if ($info[6]>0);
	$sites{$key}="@sub_sites";
	$len=$key;
	@sub_sites=();
	my %base_num=(
					$info[3] => $aa,
					$info[4] => $t,
					$info[5] => $g,
					$info[6] => $cc,
					);
	my $max_num=max(keys %base_num);
	$key_base{$key}=$base_num{$max_num};
}

my @sorted_keys=sort{$a<=>$b}keys %key_base;
my $min=$sorted_keys[0];
foreach (@sorted_keys){
	$sequence .= $key_base{$_};
}

################################################################
#	Get the location of reference
################################################################
open REF, "$ARGV[1]"  or die "$!";
my (%ref_title,%GI_seq);
$/=">";<REF>;$/="\n";
while (<REF>) {
	chomp;
	my $ref_id = $_;
	my $gi = $1 if ($ref_id =~ /gi\|(\d+)\|/);
	my $gene = (split /\s+/,$ref_id)[1];
	$ref_title{$gi}{$gene} = $ref_id;
	$/=">";
	my $ref_seq = <REF>;
	chomp $ref_seq;
	$GI_seq{$gi}{$gene} = $ref_seq;
	$/="\n";
}
close REF;

###############################################################
#	Get the 20bp from the head of CDS for next step 
###############################################################
my $hash2 = $ref_title{$GI};
my $hash3 = $GI_seq{$GI};
for $key(keys %$hash2){
	my $start_ref = (split /\s+/,$ref_title{$GI}{$key})[-2];
	my $end_ref = (split /\s+/,$ref_title{$GI}{$key})[-1];
	my $cds_len = $end_ref - $start_ref ;
	my $motify = substr($GI_seq{$GI}{$key},0,20);  # to get a motify to align with consensus sequence 

	##############################################################
	#	Get the motify location in consensu sequence
	##############################################################
#	Approximate string matching using dynamic programming
#    Find the closest match to the pattern in the text
#	This program is first condider gap than mismatch.
#	--query must shorter than --subject

	my ($pattern,$text,$gap_punishment,$mis_match);
	$gap_punishment ||= 2;
	$mis_match ||= 1;
	$text = $sequence;
	$pattern = $motify;

	my $TLEN = length $text;
	my $PLEN = length $pattern;

	# D is the Distance matrix, which shows the "edit distance" between
	# substrings of the pattern and the text.
	# It is implemented as a reference to an anonymous array.
	my $D = [  ];

	# The rows correspond to the text
	# Initialize row 0 of D.
	for (my $t=0; $t <= $TLEN ; ++$t) {
		$D->[$t][0] = 0;
	}

	# The columns correspond to the pattern
	# Initialize column 0 of D.
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
	# The text and pattern may or may not match at this character ...
					substr($text, $t-1, 1) eq substr($pattern, $p-1, 1) 
					? $D->[$t-1][$p-1]  # If they match, no increase in edit distance!
					:  $D->[$t-1][$p-1] + $mis_match,

	# gap +1 mismatch +1 
	# Second alternative
	# If the text is missing a character
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

	# Report the best match(es).
	#print "\nThe best match for the pattern $pattern\n";
	#print "has an edit distance of $bestscore\n";
	#print "and appears in the text ending at location";
	#print "s" if ( @matches > 1);
	#print " @matches\n";

	sub min3 {
		my($i, $j, $k) = @_;
		my($tmp);

		$tmp = ($i < $j ? $i : $j);
		$tmp < $k ? $tmp : $k;
	}

	# Print out
	my $p_out;my $t_out;
	foreach my $matches (@matches){
		my $plen=$PLEN;my $matches_ori=$matches;
		$p_out='';$t_out='';
		$t_out=(substr ($text, $matches-1, 1)).$t_out;
		$p_out=(substr ($pattern, $plen-1, 1)).$p_out;
		while($matches>1 && $plen>1){
			($text,$pattern,$matches,$plen,$t_out,$p_out,$D)= &add_world ($text,$pattern,$matches,$plen,$t_out,$p_out,$D);
		}
		chomp $p_out;
#	print  "pattern_location:$plen-$PLEN $p_out text_location:$matches-$matches_ori	$t_out\n";	
		my $seq_start = $matches + $min -1;
		my $terminal = $seq_start + $cds_len;
		print "$file\t$key\t$seq_start\t$terminal\n";
	}


	
sub add_world {
	my ($text_s,$pattern_s,$matches_s,$PLEN_s,$t_out_s,$p_out_s,$D_s)=@_;
	my $flag=0;
# $flag=1 <-      2 |     3 \
	if(substr ($pattern_s, $PLEN_s, 1) ne substr ($text_s, $matches_s, 1)){
		if($flag==0 && $D_s->[$matches_s][$PLEN_s-1]+$gap_punishment <= $D_s->[$matches_s-1][$PLEN_s-1]+$mis_match && $D_s->[$matches_s][$PLEN_s-1] <= $D_s->[$matches_s-1][$PLEN_s]){
			$PLEN_s=$PLEN_s-1;
			$t_out_s="_".$t_out_s;
			$p_out_s=(substr ($pattern_s, $PLEN_s-1, 1)).$p_out_s;
			$flag=1;
		}


		if($flag==0 && $D_s->[$matches_s-1][$PLEN_s] <= $D_s->[$matches_s][$PLEN_s-1] && $D->[$matches_s-1][$PLEN_s]+$gap_punishment <= $D_s->[$matches_s-1][$PLEN_s-1]+$mis_match){
			$matches_s=$matches_s-1;
			$t_out_s=(substr ($text_s, $matches_s-1, 1)).$t_out_s;
			$p_out_s="_".$p_out_s;
			$flag=2;
		}
		if($flag==0 && $D_s->[$matches_s-1][$PLEN_s-1]+$mis_match <= $D_s->[$matches_s][$PLEN_s-1]+$gap_punishment && $D_s->[$matches_s-1][$PLEN_s-1]+$mis_match <= $D_s->[$matches_s-1][$PLEN_s]+$gap_punishment){
			$matches_s=$matches_s-1;
			$PLEN_s=$PLEN_s-1;
			$t_out_s=(substr ($text_s, $matches_s-1, 1)).$t_out_s;
			$p_out_s=(substr ($pattern_s, $PLEN_s-1, 1)).$p_out_s;
			$flag=3;
		}
	}
	if(substr ($pattern_s, $PLEN_s, 1) eq substr ($text_s, $matches_s, 1)){
		if($flag==0 && $D_s->[$matches_s][$PLEN_s-1]+$gap_punishment <= $D_s->[$matches_s-1][$PLEN_s-1] && $D_s->[$matches_s][$PLEN_s-1] <= $D_s->[$matches_s-1][$PLEN_s]){
			$PLEN_s=$PLEN_s-1;
			$t_out_s="_".$t_out_s;
			$p_out_s=(substr ($pattern_s, $PLEN_s-1, 1)).$p_out_s;
			$flag=1;
		}
		if($flag==0 && $D_s->[$matches_s-1][$PLEN_s] <= $D_s->[$matches_s][$PLEN_s-1] && $D->[$matches_s-1][$PLEN_s]+$gap_punishment <= $D_s->[$matches_s-1][$PLEN_s-1]){
			$matches_s=$matches_s-1;
			$t_out_s=(substr ($text_s, $matches_s-1, 1)).$t_out_s;
			$p_out_s="_".$p_out_s;
			$flag=2;
		}
		if($flag==0 && $D_s->[$matches_s-1][$PLEN_s-1] <= $D_s->[$matches_s][$PLEN_s-1]+$gap_punishment && $D_s->[$matches_s-1][$PLEN_s-1] <= $D_s->[$matches_s-1][$PLEN_s]+$gap_punishment){
			$matches_s=$matches_s-1;
			$PLEN_s=$PLEN_s-1;
			$t_out_s=(substr ($text_s, $matches_s-1, 1)).$t_out_s;
			$p_out_s=(substr ($pattern_s, $PLEN_s-1, 1)).$p_out_s;
			$flag=3;
		}
	}
	return ($text_s,$pattern_s,$matches_s,$PLEN_s,$t_out_s,$p_out_s,$D_s) if($flag!=0);
}


}







