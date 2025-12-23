#!/usr/bin/perl

=head1 Name
	
	primer_like.pl - use a block to hit conserved regions of the sequences just
	like use a series of primers to extract expected DNA sequences.

=head1 Description
	This script has two parts, the core program is writed by Hailin Pan who is 
	doing research in BGI, and rest part is writed by Chentao Yang.
	Approximate string matching using dynamic programming
	Find the closest match to the pattern in the seq
	This program is first condider gap than mismatch.
	--query must shorter than --subject

=head1 Version
	
	Author: Chentao Yang, yangchentao@genomics.cn
	Version: 1.0,  Date: 2015-02-1

=head1 Attention

	This program must be used with "string_match.pl", you can copy that into 
	the same directory which this program is in.

=head1 Usage
	
	perl $0  [option]
	--primer <str> query string
	--fa <str>	sequence fasta format
	--gap_punish <num> gap punish	default value 2		
	--mismatch_punish <num> mismatch punish	default value 1
	--cutoff_mis <num> mismatch allowable 	default value 4
	--cutoff_gap <num> gap allowable	default value 2
	--log <str> output file
	--out <str>	output the sequences

=head1 Exmple

	perl ./primer_like.pl --primer primer.txt --fa seq.txt -o output_file &
	perl ./primer_like.pl --primer primer.txt --fa seq.txt --gap_punish 4 --mismatch_punish 2 -cutoff_mis 4 -cutoff_gap 2 -log extract.seq.txt.log -o  output_file 
	primer.txt may be like this: EIQADEVRL 
	seq.txt may be like this: 
	>id1
	SVLQDRSMPHQEILAADEVLQESEMRQQDMISHDE
=cut
use Cwd;
use strict;
use warnings;
use Getopt::Long;

my ($primer,$fa,$cutoff_mis,$cutoff_gap,$output,$log,$help);
GetOptions(
	"primer:s"=>\$primer,
	"fa:s"=>\$fa,
	"cutoff_mis:i"=>\$cutoff_mis,
	"cutoff_gap:i"=>\$cutoff_gap,
	"log:i"=>\$log,
	"out:s"=>\$output,
	"help!"=>\$help
);

$cutoff_mis ||=4;
$cutoff_gap ||=2;
$log ||="extract.$fa.log";
die `pod2text $0` if (!($primer && $fa)|| $help);
#Attention ,this commond will remove the file "$output",make sure that you do not have the same name file;
`rm $log if [ -f $log ]`;
open QUE,"$primer";
open OBJ,"$fa";
open LOG, ">>$log";
my $workdir = getcwd();
my $pattern = <QUE>;
chomp $pattern;
#Loop each sequence for check
my %hash;
$/=">";<OBJ>;$/="\n";
while (my $id = <OBJ>) {
	chomp $id;
	$/=">";
	my $seq = <OBJ>;
	chomp $seq;
	$seq =~ s/\n//g;
	$hash{$id} = $seq;
	print LOG "//$id\n";
	`perl $workdir/string_match.pl --query "$pattern" --subject "$seq" >>$log`;
	$/="\n";
	
}
close QUE;
close OBJ;
close LOG;

###check out result
open IN, "$log";
open OUT, ">$output";
$/="//";<IN>;$/="\n";
while (my $id = <IN>) {
	chomp $id;
	$/="//";
	my $len = <IN>;
	chomp $len;
	my @arr = split /\n/,$len;
	my $mark = 0;
	foreach my $result(@arr) {
		my $p_out = (split /\s+/,$result)[1];
		my $t_out = (split /\s+/,$result)[3];
		my $len_out = length($t_out);
		my $mismatch_count = 0;my $gap_count = 0;
		for (my $x = 0,my $y = 0;$x<=$len_out-1,$y<=$len_out-1;$x++,$y++) {
			my $tmp_x = substr($p_out,$x,1);
			my $tmp_y = substr($t_out,$y,1);
			$mismatch_count++ if (($tmp_x ne "_") && ($tmp_y ne "_") && ($tmp_x ne $tmp_y)) ;
			$gap_count++ if ($tmp_y eq "_" || $tmp_x eq "_");
		}
		print "$mismatch_count,$gap_count\n";
		if (($mismatch_count <= $cutoff_mis) && ($gap_count <= $cutoff_gap)) {
			$mark++;
			$mismatch_count = 0;$gap_count = 0;
		}
	}
	if ($mark>0){ print OUT ">$id\n$hash{$id}\n" ;$mark = 0 ;}

	$/="\n";
}
close IN;
