#!/usr/bin/perl -w
use strict;
use List::Util qw(max min);

die "Usage: perl $0 <sort.gff> <repeat type (DNA|LTR|SINE|LINE...)>\n" if @ARGV <2;

my $file = shift;
my $repeat = shift;
open IN,"<$file" or die "$!\n";
open OUT,">$file.combine";

my (@types,@starts,@ends);
my $str = <IN>;
chomp $str ;
my @a = split /\s+/,$str;
my $tmp_name = $a[0];
push @types,$a[1];
my $tmp_start = $a[2];
my $tmp_end = $a[3];
push @starts,$tmp_start;
push @ends,$tmp_end;

while (<IN>) {
	next if ($_ eq $str);
	next if (/#/);
	chomp;
	my @t = split /\s+/;
	my $name = $t[0];
	my $type = $t[1];
	my $start = $t[2];
	my $end = $t[3];
	die "input file is not sorted, which means end position is less than start position!" 
	unless ($end > $start) ;
	if($start > max(@ends)){
		my %count;
		my @uniq_types = grep { ++$count{ $_ } < 2; } @types;
		my $num = @types;
		if (@uniq_types == 1){
			my $st = min(@starts);
			my $en = max(@ends);
			print OUT "$tmp_name\t$types[0]\t$st\t$en";
			$num == 1 ? print OUT "\tINNOCENT\t-\n" : print OUT "\tCombined\t-\n";
			@types = ();
			push @types,$type;
			@starts = ();
			push @starts,$start;
			@ends = ();
			push @ends,$end;
			$tmp_name = $name;

		}else {
			my $ty = join("|",@types);
			my $st = min(@starts);
			my $en = max(@ends);
			print OUT "$tmp_name\t$repeat\t$st\t$en\tCombined(Complicated)\t$ty\n";
			@types = ();
			push @types,$type;
			@starts = ();
			push @starts,$start;
			@ends = ();
			push @ends,$end;
			$tmp_name = $name;
		}
	
	}elsif( $name ne $tmp_name) {
		my %count;
		my @uniq_types = grep { ++$count{ $_ } < 2; } @types;
		my $num = @types;
		if (@uniq_types == 1){
			my $st = min(@starts);
			my $en = max(@ends);
			print OUT "$tmp_name\t$types[0]\t$st\t$en";
			$num == 1 ? print OUT "\tINNOCENT\t-\n" : print OUT "\tCombined\t-\n";
			@types = ();
			push @types,$type;
			@starts = ();
			push @starts,$start;
			@ends = ();
			push @ends,$end;
			$tmp_name = $name;

		}else {
			my $ty = join("|",@types);
			my $st = min(@starts);
			my $en = max(@ends);
			print OUT  "$tmp_name\t$repeat\t$st\t$en\tCombined(Complicated)\t$ty\n";
			@types = ();
			push @types,$type;
			@starts = ();
			push @starts,$start;
			@ends = ();
			push @ends,$end;
			$tmp_name = $name;
		}
	}elsif($start < max(@ends)){
		push @types,$type;
		push @starts,$start;
		push @ends,$end;

	}

}

close IN;
