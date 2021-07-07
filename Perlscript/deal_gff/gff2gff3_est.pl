#!/usr/bin/perl
die"usage: perl gff2gff3_est.pl <gff file> <source name>[species_name|gene_prediction_tools_name|...]\n"unless @ARGV >= 1;
my $files = @ARGV;
my ($gff)= shift;

my ($source,$feature);
if ($files > 1) {
	$source= shift;
	$feature="$source\_match";
}


open(IN, "<$gff")||die"can not open $gff\n";
while(<IN>) {
  my @t=split/\s+/,$_;
  my $id=$t[8];
  my @p=split/\../,$t[9];
  $t[8]="ID=match.$id;Target=$id $p[0] $p[1]";
  if(defined $source) {
  	print"$t[0]\t$source\t$feature\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\t$t[8]\n";
  }else {
  	$feature = "$t[1]\_match";
  	print"$t[0]\t$t[1]\t$feature\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\t$t[8]\n";
  }
  
}

close IN;
