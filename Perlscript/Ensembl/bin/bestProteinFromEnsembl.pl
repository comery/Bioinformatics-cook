#!/usr/bin/perl -w
use strict;
die "Usage: <pep.all.fa.gz (Ensembl format)> <gff.orf> <priority[1:intact ORF; 2:length]>\n" unless @ARGV == 3;

die unless $ARGV[2] == 1 || $ARGV[2] == 2;

my %intact;
if ($ARGV[1] =~ /\.gz$/) {
	open IN, "gunzip -c $ARGV[1] | ";
} else {
	open IN, $ARGV[1];
}
while (<IN>) {
	my @info = split /\s+/;	
	$intact{$info[0]} += $info[7] + $info[8];# if $info[7] && $info[8];
}
close IN;

my %pep;
if ($ARGV[0] =~ /\.gz$/) {
	open IN, "gunzip -c $ARGV[0] | ";
} else {
	open IN, $ARGV[0];
}
$/ = ">";
<IN>;
while (<IN>) {
	/(.+)\n/;
	my $head = $1;
	my $protein_id = (split /\s+/, $head)[0];
	next unless defined($intact{$protein_id});
	die unless $head =~ /gene:(\S+?)\s+/;
	my $gene_id = $1;
	s/.+\n//;
	s/>|\s+//g;
	my $len = length($_);
	push @{$pep{$gene_id}}, [$protein_id, $len, $intact{$protein_id}];
} 
$/ = "\n";
close IN;


foreach my $gene_id (keys %pep) {
	@{$pep{$gene_id}} = sort {$b->[-1] <=> $a->[-1] or $b->[-2] <=> $a->[-2]} @{$pep{$gene_id}} if $ARGV[2] == 1;	
	@{$pep{$gene_id}} = sort {$b->[-2] <=> $a->[-2] or $b->[-1] <=> $a->[-1]} @{$pep{$gene_id}} if $ARGV[2] == 2;	
	foreach my $p (@{$pep{$gene_id}}) {
		my ($protein_id, $len, $intact_num) = @$p;
		print "$protein_id\t$gene_id\t$len\t$intact_num\n";
		last;
	}
}
