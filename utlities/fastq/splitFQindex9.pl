#!/usr/bin/perl
# Split fq file according to index(head 8bp) and the GATC(AATTC)

use strict;
use List::Util qw(max min);
use PerlIO::gzip;

die "Usage: perl $0 <index.lst> <lane.lst> <fq1> <fq2> <outdir> <mis num> <column SID> <column Enzyme>\n" if(@ARGV != 8);
my $findex = shift;
my $flane = shift;
my $ffq1 = shift;
my $ffq2 = shift;
my $outdir = shift;
my $misnum = shift;
my $sampleID = shift;
my $colEnzyme = shift;

`mkdir -p $outdir` unless(-d $outdir);

my $prefix = $1 if($ffq1 =~ /.*(\d{6}_\w\d{3}_\w+_L\d+_[^_]+)_1\./);
if(!$prefix){
	$prefix = "none";
}

my (@maxindex,$maxindex,,$minindex,%index);
open IN,$findex or die $!;
while(<IN>){
	chomp;split;
	$index{$_[0]} = $_[1];
	push @maxindex,length($_[1]);
}
close IN;
$maxindex = max @maxindex;
$minindex = min @maxindex;
print STDERR "Max and Min index length : $maxindex, $minindex\n";

# hash table
my @base = ('A','C','G','T','N');
my ($lane,%fh,%real,%hash1,%hash2,%has1,%has2,%hass1,%hass2);
open IN,$flane or die $!;
while(<IN>){
	chomp;split;
	if($_[10] =~ /^\w/){
		$lane = $_[10];
		open OM,">:gzip","$outdir/$lane.mis.reads.gz";
	}
	`mkdir -p $outdir/$_[$sampleID - 1]`;
	if(!$index{$_[$colEnzyme - 1]}){ die "Error: lack info for enzyme $_[$colEnzyme-1]!\n"; }
	my $indexlen = length($index{$_[$colEnzyme-1]});
	my $seq = $index{$_[$colEnzyme-1]}.$index{enzyme};
	$seq = $1 if($seq =~ /^(\w{8})/);
	$real{$seq}{seq} = $_[$sampleID - 1];
	$real{$seq}{len} = $indexlen;
	if($misnum >= 1){
		hashTable(\%hash1,$seq,1,$_[$sampleID - 1],$indexlen);
		foreach(keys %hash1){
			if(!$has1{$_}){
				$has1{$_} = $hash1{$_}{seq};
				$hass1{$_}{$hash1{$_}{seq}} = 1;
			}elsif(!$hass1{$_}{$hash1{$_}{seq}}){
				$has1{$_}.= ":".$hash1{$_}{seq};
				$hass1{$_}{$hash1{$_}{seq}} = 1;
			}
		}
	}

	if($misnum >= 2){
		hashTable(\%hash2,$seq,2,$_[$sampleID - 1],$indexlen);
		foreach(keys %hash2){
			if(!$has2{$_}){
				$has2{$_} = $hash2{$_}{seq};
				$hass2{$_}{$hash2{$_}{seq}} = 1;
			}elsif(!$hass2{$_}{$hash2{$_}{seq}}){
				$has2{$_}.= ":".$hash2{$_}{seq};
				$hass2{$_}{$hash2{$_}{seq}} = 1;
			}
		}
	}

	open $fh{$_[$sampleID - 1]}{1},">:gzip","$outdir/$_[$sampleID - 1]/$_[$sampleID - 1]-$prefix\_1.fq.gz";
	open $fh{$_[$sampleID - 1]}{2},">:gzip","$outdir/$_[$sampleID - 1]/$_[$sampleID - 1]-$prefix\_2.fq.gz";
}
close IN;

if($misnum >= 1){
	foreach(keys %has1){
		if($has1{$_} =~ /:/){
			delete $hash1{$_};
		}
	}
}

if($misnum >= 2){
	foreach(keys %has2){
		if($has2{$_} =~ /:/){
			delete $hash2{$_};
		}
	}
}

open I1,"<:gzip",$ffq1 or die $!;
open I2,"<:gzip",$ffq2 or die $!;

my ($line1,$line2,$hseq,@fq1,@fq2);
my $count = 0;
while($line1 = <I1> and $line2 = <I2>){
	$count++;
	if($count <= 4){
		push @fq1,$line1; push @fq2,$line2;
	}
	if($count == 2){
		$hseq = $1 if($line1 =~ /^(\w{8})/);
	}
	if($count == 4){
		$count = 0;
		my ($restN,$N) = (0,0);
		if($real{$hseq}{seq}){ # no mismatch
			$N = $real{$hseq}{len};
			$restN = $maxindex - $N;
			if($fq1[1] =~ /^[ACGTN]{$N}$index{enzyme}/){
				($fq1[1],$fq1[3]) = trimRead($fq1[1],$fq1[3],$N,$restN);
				my ($fh1,$fh2) = ($fh{$real{$hseq}{seq}}{1},$fh{$real{$hseq}{seq}}{2});
				print $fh1 @fq1; print $fh2 @fq2;
			}else{
				print OM "#$hseq\n"; print OM @fq1,@fq2;
			}
		}elsif($hash1{$hseq}{seq}){ # one mismatch
			$N = $hash1{$hseq}{len};
			$restN = $maxindex - $N;
			if($fq1[1] =~ /^[ACGTN]{$N}$index{enzyme}/){
				($fq1[1],$fq1[3]) = trimRead($fq1[1],$fq1[3],$N,$restN);
				my ($fh1,$fh2) = ($fh{$hash1{$hseq}{seq}}{1},$fh{$hash1{$hseq}{seq}}{2});
				print $fh1 @fq1; print $fh2 @fq2;
			}else{
				print OM "#$hseq\n"; print OM @fq1,@fq2;
			}
		}elsif($hash2{$hseq}{seq}){ # two mismatch
			$N = $hash2{$hseq}{len};
			$restN = $maxindex - $N;
			if($fq1[1] =~ /^[ACGTN]{$N}$index{enzyme}/){
				($fq1[1],$fq1[3]) = trimRead($fq1[1],$fq1[3],$N,$restN);
				my ($fh1,$fh2) = ($fh{$hash2{$hseq}{seq}}{1},$fh{$hash2{$hseq}{seq}}{2});
				print $fh1 @fq1; print $fh2 @fq2;
			}else{
				print OM "#$hseq\n"; print OM @fq1,@fq2;
			}
		}else{
			print OM "#$hseq\n"; print OM @fq1,@fq2;
		}
		undef @fq1; undef @fq2;
	}
}
close I1; close I2;
print STDERR "Done.\n";

#===========#
sub trimRead{
#===========#
	my ($read,$qv,$Nh,$Nt) = @_;
	$read =~ s/^\w{$Nh}//;
	$qv =~ s/^.{$Nh}//;
	if($Nt > 0){
		$read =~ s/\w{$Nt}$//;
		$qv =~ s/.{$Nt}$//;
	}
	return ($read,$qv);
}

#============#
sub hashTable{
#============#
	my ($p,$seq,$mis,$v,$indexlen) = @_;
	my ($i,$j) = (0,0);
	my $len = length($seq);
	for($i = 0; $i < $len; ++$i){
		for(my $ii = 0; $ii < @base; ++$ii){
			my $seq1 = $seq;
			substr($seq1,$i,1) = $base[$ii];
			if($mis == 1){
				$p->{$seq1}{seq} = $v;
				$p->{$seq1}{len} = $indexlen;
			}
			if($mis == 2){
				for($j = $i+1; $j < $len; ++$j){
					for(my $ii = 0; $ii < @base; ++$ii){
						my $seq2 = $seq1;
						substr($seq2,$j,1) = $base[$ii];
						$p->{$seq2}{seq} = $v;
						$p->{$seq2}{len} = $indexlen;
					}
				}
			}
		}
	}
}

