#!/usr/bin/perl -w
use strict;
use lib '/share/project002/liqiye/bin/module/personal';
use GFF;
use Fasta;
die "Usage: <gff> <ref.fa>\n" unless @ARGV == 2;
################################################################
# detect orf without extending. Author: liqiye@genomics.org.cn #
################################################################

my %gene_pos;
getChrGene($ARGV[0], \%gene_pos);

my %starCodon = ("ATG"=>1);
my %stopCodon = ("TAA"=>1, "TGA"=>1, "TAG"=>1);
if ($ARGV[1] =~ /\.gz$/) {
	open IN, "gunzip -c $ARGV[1] | ";
} else {
	open IN, $ARGV[1];
}
$/ = ">";
<IN>;
while (<IN>) {
	/(.+)\n/;
	my $chr = (split /\s+/, $1)[0];
	next unless $gene_pos{$chr};
	s/.+\n//;
	s/\s+|>//g;
	my $len = length($_);
	foreach my $a_p (@{$gene_pos{$chr}}) {
		my ($gene, $bg, $ed, $strand) = @$a_p;
		my ($star, $stop);
		if ($strand eq "+") {
			$star = substr($_, $bg-1, 3);
			$stop = substr($_, $ed-3, 3);
			unless ($stopCodon{$stop}) {
				unless ($ed == $len) {
					$stop = substr($_, $ed, 3);
					$ed += 3;
				}
			}
		} else {
			$star = substr($_, $ed-3, 3);
			$star = reverse_complement($star);
			$stop = substr($_, $bg-1, 3);
			$stop = reverse_complement($stop);
			unless ($stopCodon{$stop}) {
				unless ($bg == 1) {
					$stop = substr($_, $bg-4, 3);
					$stop = reverse_complement($stop);
					$bg -= 3;
				}
			}
		}
		$star = uc($star);
		$stop = uc($stop);
		my ($is_star, $is_stop);
		$is_star = $starCodon{$star} ? 1 : 0;
		$is_stop = $stopCodon{$stop} ? 1 : 0;
		print "$gene\t$chr\t$strand\t$bg\t$ed\t$star\t$stop\t$is_star\t$is_stop\n";
	}
}
$/ = "\n";
close IN;
