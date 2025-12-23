#!/usr/bin/perl
  die"usage: perl gff2gff3_est.pl <gff file>\n"unless @ARGV==1;
  my ($gff)=@ARGV;
  my $source="EST";
  my $feature="EST_match";
  open(IN, "<$gff")||die"can not open $gff\n";
  while(<IN>)
  {
      my @t=split/\s+/,$_;
      $t[1]=$source;
      $t[2]=$feature;
      my $id=$t[8];
      my @p=split/\../,$t[9];
      $t[8]="ID=match.$id;Target=$id $p[0] $p[1]";
      print"$t[0]\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\t$t[8]\n";
}
close IN;