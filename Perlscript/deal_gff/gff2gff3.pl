#!/usr/bin/perl -w

use strict;
#use DataUtils;
#use CombUtils;
use FileHandle;

#######################################################
## to use EVAL scripts, make a gtf file; convert gff 2  gff3
#######################################################
if (@ARGV < 1){
	print "Usage: perl $0 <input.gff>";
	exit;
}


gff2Gff3();

sub gff2Gff3 {
  my %genes=loadGff3P();

  my $gcnt=0;
  my %nLabel;
  print "##gff-version 3\n";
  foreach my $seqId (keys %genes) {
    my $gLst=$genes{$seqId};
    my %ghsh=%$gLst;
    my @nGenes=();
    my @pGenes=();
    my %oldLabel;
    foreach my $gid (keys %ghsh) {
      my $gene = $ghsh{$gid};
      my @gc=split(/:/,$gene);
      #print "@gc\n";
      my $strnd = shift @gc;
      my @sc = sort { (split(/ /,$a))[0] <=> (split(/ /,$b))[0] } @gc;
      my ($end5,$ignore,$ignore2,$ignore3,$pid) = split(/ /,$sc[0]);
      my ($ignoreA,$end3,$ignore4,$ignore5,$pid3) = split(/ /,$sc[$#sc]);
      my @mrna=("$seqId",$pid,"mRNA",$end5,$end3,".",$strnd,".","ID=$gid.gene;Name=$gid.gene");
      prnLine(\@mrna);
      my $exnum=1;
      for(my $iter = 0; $iter <= $#gc; $iter += 1) { 
	my ($end5,$end3,$phase,$type,$pid) = split(/ /,$gc[$iter]);
	my @cds=("$seqId",$pid,"CDS",$end5,$end3,".",$strnd,$phase,"ID=$gid.cds.$exnum;Parent=$gid.gene;Name=$gid.gene;Note=$type");
	prnLine(\@cds);
	$exnum++;
      }
    }
  }
}

sub prnLine {
  my ($lRef) = @_;
  my @line=@$lRef;
  print "$line[0]";
  for(my $iter = 1; $iter <= $#line; $iter++) {
    print "\t$line[$iter]";
  }
  print "\n";
}


sub loadGff3P {
  my %genes=();
  my $strndIdx = 6;
  my ($idx1,$idx2) = (3,4);
  open IN, $ARGV[0];
  while (my $line=<IN>) {
    chomp($line);
    if ($line =~ /\#/ || !$line ) {
      next;	# skip coment line
    }
    my @vals = split(/\t/,$line);
    my $seqname = $vals[0];
    my $strnd = $vals[$strndIdx];
    my $gid = $vals[$#vals];
    if( !$seqname ) {
	print "huh? [$line]\n";
    }
    if( !$genes{$seqname} || !$genes{$seqname}{$gid} ) {
      $genes{$seqname}{$gid} = "$strnd:";
      if( !$strnd ) {
	print "error $line\n";
	exit(0);
      }
    }
    $genes{$seqname}{$gid} .= "$vals[$idx1] $vals[$idx2] $vals[7] $vals[2] $vals[1]:";
  }
  return %genes;
  close IN;
}

