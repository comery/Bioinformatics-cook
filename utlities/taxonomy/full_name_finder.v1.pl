#! usr/bin/perl -w
use strict;
use FindBin qw($Bin $Script);
die "perl <species name in the format of genus_species>" unless (@ARGV==1);
open NOD, "$Bin/taxdmp/nodes.dmp" || die $!;
my (%up,%lev);
while (<NOD>){
	chomp;
	my @a = split (/\|/,$_,4);
	$a[0]=~s/\s+//g;
	$a[1]=~s/\s+//g;
	$a[2]=~s/\s+//g;
	$up{$a[0]}=$a[1];
	$lev{$a[0]}=$a[2];
}
close NOD;
open NAM, "$Bin/taxdmp/names.dmp" || die $!;
my %name;
my %hash;
while (<NAM>){
	chomp;
	my @in=split (/\|/,$_,3);
	$in[0]=~s/\s+//g;
	$hash{$in[1]}=$in[0];
	next unless (/scientific name/);
	my @a = split (/\|/,$_,3);
	$a[0]=~s/\s+//g;
	$a[1]=~s/^\s+//;
	$a[1]=~s/\s+$//;
	$name{$a[0]}=$a[1];
}
close NAM;
open IN, "$ARGV[0]" || die $!;
open OUT, ">$ARGV[0]\.full" || die $!;
while(<IN>){
	chomp;
	my @genful= split /\_/;
	my $gen="$genful[0]"." ";
#	my $genus="$genful[0]" if ($genful[1]=~/sp./);
#	print "don't have $_ in database, it could be error\n" unless (exists $rname{$_});
        my $neo;
	my $jd="genus";
	my $jug=0;
	NNN:for my $key (keys %name){
		if ($name{$key} eq $genful[0] and $lev{$key} eq $jd){
			$neo=$key;
			$jug=1;
			last NNN;
		}
	}
	if ($jug==0){
		TTT:for my $key (sort keys %hash){
			if ($key=~/$gen/){
				$neo=$hash{$key};
				last TTT;
			}
		}
	}
	my @all;
#	print "!!!$neo\t$hash{$neo}\n";
	unless ($neo){
		print "error in $_, $_ can't be found in the database\n";
		next;
	}
	while(1){
		if($lev{$neo} eq "kingdom"){
                        $all[0] = $neo;
                }elsif($lev{$neo} eq "phylum"){
                        $all[1] = $neo;
                }elsif($lev{$neo} eq "class"){
                        $all[2] = $neo;
                }elsif($lev{$neo} eq "order"){
                        $all[3] = $neo;
                }elsif($lev{$neo} eq "suborder"){
                        $all[4] = $neo;
                }elsif($lev{$neo} eq "superfamily"){
                        $all[5] = $neo;
                }elsif($lev{$neo} eq "family"){
                        $all[6] = $neo;
                }elsif($lev{$neo} eq  "subfamily"){
                        $all[7] = $neo;
                }elsif($lev{$neo} eq "genus"){
                        $all[8] = $neo;
                }
		my $int=$up{$neo};
		$neo=$int;
		last if $neo==1;
	}
	print OUT "$_\t";
	for my $i (0..8){
		if (defined $all[$i]){
			print OUT "$name{$all[$i]}\t"
		}else{
			print OUT "NA\t"
		}
	}
	print OUT "\n";
}
