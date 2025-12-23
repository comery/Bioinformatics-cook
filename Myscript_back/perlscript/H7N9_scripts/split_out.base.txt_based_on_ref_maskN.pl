#!/usr/bin/perl -w
use strict;
use List::Util qw(max);
#######################################
#
# Set up usage statement
#
#######################################
my $usage = "Warnning:\n\t Usage: perl  split_out.base.txt_based_on_ref_maskN.pl <out.base.txt> <cluster_ref.fasta.cor.maskN.fa> <outdir>\n";

my $scripthelp = "\tsplit_out.base.txt_based_on_ref_maskN.pl - split out.base.txt file into different segment based on best reference\n";

#######################################
#
# Initialize variables
#
#######################################
my $outbase;
my $ref;
my $outdir;
my $argNum = 3;

#######################################
#
# Test for commandline arguments
#
#######################################

if (! $ARGV[0] ) {
	print $usage;
	print $scripthelp;
	exit -1;
} elsif (scalar @ARGV != $argNum) {
	print $usage;
	print $scripthelp;
	exit -1;
} 

$outbase = $ARGV[0];
if (! -f $outbase) {
	print "Unable to locate input fasta file: $outbase.\n";
	exit;
}
$ref = $ARGV[1];
$ARGV[2] ? $outdir = $ARGV[2] :$outdir = "outdir";
`mkdir $outdir` if (! -e $outdir);
chdir "$outdir";
`mkdir seg1 seg2 seg3 seg4 seg5 seg6 seg7 seg8`;
chdir "../";
#########################################
#
# Open the file
#
########################################
open REF,"<$ref";
my (%hash,$ref_id,$id,$ref_len,$coverage,$true_len,$GI,$seg,$ref_seq);
$/=">";<REF>;$/="\n";
while ($ref_id = <REF>) {  ## make a hash to record the N number and gi number
	chomp $ref_id;
	my @aa;
	if ($ref_id =~ /:\_/) {
		@aa = split /\_/,$ref_id;  ## sometimes you're not sure what the pattenn is ,it will help you.
	}else{
		@aa = split /\s+/,$ref_id;
	}
	$id = $aa[0];
	$ref_len = $1 if ($aa[1] =~ /length(\d+)/);
	$coverage = $aa[-1];
	$true_len = $ref_len * $coverage;
	$GI = $1 if ($id =~ m/gi\|(\d+)\|gb/);
	$seg = $2 if ($id  =~ m/Influenza-(A-)?(\w{3}\d)\_?/);
	$seg = lc($seg);
	$/=">";
	$ref_seq = <REF>;
	chomp $ref_seq;
#	if ($ref_seq =~ /N/) {
#		$N_num = ($ref_seq =~ s/N/N/g);
#	}else{
#		$N_num = 0;
#	}
	$hash{$seg}{$true_len} = $GI;
	$/="\n";
}

my %hash3;
open ID,">id.info.list";
foreach my $key1 (keys %hash){
	my $hash2 = $hash{$key1};
	my $max_len = max(keys %$hash2);
	$hash3{$key1} = $hash{$key1}{$max_len};
	print ID "$key1\t$hash3{$key1}\t$max_len\n" ; ## print the best reference's gi number and it's real coverage
}
close ID;
close REF;
open BASE,"<$outbase";
open OUT1,">$outdir/seg1/seg1-out.base.txt";
open OUT2,">$outdir/seg2/seg2-out.base.txt";
open OUT3,">$outdir/seg3/seg3-out.base.txt";
open OUT4,">$outdir/seg4/seg4-out.base.txt";
open OUT5,">$outdir/seg5/seg5-out.base.txt";
open OUT6,">$outdir/seg6/seg6-out.base.txt";
open OUT7,">$outdir/seg7/seg7-out.base.txt";
open OUT8,">$outdir/seg8/seg8-out.base.txt";

while (<BASE>) {	## sleect the different segment's data by gi number which is best (N number is the least!)
	next if ($_ =~ /Scaffold/);
	print OUT1 "$_" if ($_ =~ /$hash3{seg1}/) ;
	print OUT2 "$_" if ($_ =~ /$hash3{seg2}/) ;
	print OUT3 "$_" if ($_ =~ /$hash3{seg3}/) ;
	print OUT4 "$_" if ($_ =~ /$hash3{seg4}/) ;
	print OUT5 "$_" if ($_ =~ /$hash3{seg5}/) ;
	print OUT6 "$_" if ($_ =~ /$hash3{seg6}/) ;
	print OUT7 "$_" if ($_ =~ /$hash3{seg7}/) ;
	print OUT8 "$_" if ($_ =~ /$hash3{seg8}/) ;
}

close BASE;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
close OUT7;
close OUT8;
