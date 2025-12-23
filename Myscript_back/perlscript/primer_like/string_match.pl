#!/usr/bin/perl

=head1 Name
	
	string_match.pl - find the best match of two strings

=head1 Description

	Approximate string matching using dynamic programming
	Find the closest match to the pattern in the text
	This program is first condider gap than mismatch.
	--query must shorter than --subject

=head1 Version
	
	Author: HailinPan, panhailin@genomics.org.cn
	Version: 1.0,  Date: 2011-02-26

=head1 Usage
	
	perl $0  [option]
	--query <str> query string
	--subject <str>	subject string
	--gap_punish <num> gap punish
	--mismatch_punish <num> mismatch punish

=head1 Exmple

	perl ./string_match.pl --query EIQADEVRL --subject SVLQDRSMPHQEILAADEVLQESEMRQQDMISHDE > output_file &
	perl ./string_match.pl --query EIQADEVRL --subject SVLQDRSMPHQEILAADEVLQESEMRQQDMISHDE --gap_punish 4 --mismatch_punish 2 > output_file &

=cut

use strict;
use warnings;
use Getopt::Long;

my ($pattern,$text,$gap_punishment,$mis_match);
my $help;
GetOptions(
	"query:s"=>\$pattern,
	"subject:s"=>\$text,
	"gap_punish:i"=>\$gap_punishment,
	"mismatch_punish:i"=>\$mis_match,
	"help!"=>\$help
);
$gap_punishment ||= 2;
$mis_match ||= 1;
die `pod2text $0` if (!($pattern && $text)|| $help);
#$pattern = 'EIQADEVRL';
#print "PATTERN:\n$pattern\n";
#$text = 'SVLQDRSMPHQEILAADEVLQESEMRQQDMISHDE';
#print "TEXT:\n$text\n";
#print "gap_punish:$gap_punishment\n";
#print "mismatch_punish:$mis_match\n";

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
#	$D->[0][$p] = 0;
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

# Print D, the resulting edit distance array
#for (my $p=0; $p<=$PLEN; $p++){
#	print "        "if($p==0);
#	print substr($pattern,$p-1,1), "   " if($p>0);
#}
#print "\n";
#for (my $t=0; $t <= $TLEN ; ++$t) {
#	print substr($text,$t-1,1), " " if($t>0);
#	print "  "if($t==0);
#	for (my $p=0; $p <= $PLEN ; ++$p) {
#		printf "%3.0f",$D->[$t][$p]; print " ";
#	}
#	print "\n";
#}

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
	print  "pattern_location:$plen-$PLEN $p_out text_location:$matches-$matches_ori	$t_out\n";	
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










