#! /usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $usage = "Usage: $0 -f <fasta file> -o <relaxed Phylip file> \n";
$usage .= "            -f   <input fasta file>\n";
$usage .= "            -o   <relaxed Phylip file>\n";
$usage .= "            -h   Displays this help and exit\n\n";

my $fasta;
my $help;
my $outfile;

GetOptions(
           'f=s'	  => \$fasta,
           'o=s'         => \$outfile,
           'h'        => \$help,
          );

if( $help ){
    print $usage;
    exit 0;
}

unless( $fasta && -f $fasta && $outfile ){
    print STDERR $usage;
    exit 1;
}

#########################################################################
### i. Open Fasta
### ii. Store sequence (count lenght and number)
### iii. Return sequence correctly formated (one line per seq)
### 	The header should contain a space chataracter and then number of
###	sequence as well as the length of the alignement

#Change line delimiter :
$/='>' ;

open ( my $IN , "<$fasta" ) or die "Unable to open $fasta !\n " ;

my %Alignements ;
my $biggestSequence = 0 ;
while (my $i = <$IN> ){
	my @words = split "\n", $i ;

	# process header
	my $header = shift @words ;
	$header =~ s/\n//g ;
	$header =~ s/>// ;
	
	# process sequence
	my $seq = join '', @words ;
	$seq =~ s/\n//g ;
	$seq =~ s/>//g ;
	my $length = length $seq ;
	if( $length > $biggestSequence ){ $biggestSequence = $length ; }
	print "$header : $length\n" ;
	
	# store information
	$Alignements{$header}{'Length'} = $length ;
	$Alignements{$header}{'Sequence'} = $seq ;

}

close $IN ;

# Check sequence, determine number of sequence :
my $numberOfSequence = 0 ;
foreach my $i ( keys %Alignements ){
	unless( $Alignements{$i}{'Length'} == $biggestSequence ){ 
		print "$i doesn't have the same length than others ($Alignements{$i}{'Length'}/$biggestSequence)!\n " ; 
		next ;
	}
	$numberOfSequence++ ;
} 

open OUT , ">$outfile" or die "Unable to open $outfile !\n" ;
print OUT " $numberOfSequence $biggestSequence\n" ;
foreach my $i ( keys %Alignements ){
	if( $Alignements{$i}{'Length'} == $biggestSequence ){ print OUT "$i\t$Alignements{$i}{'Sequence'}\n" ;}
}
close OUT ;
