#!/usr/bin/perl

=head1 Name

getGene.pl  --  get gene elements from sequences according to coordinates

=head1 Description

This program is used to get out the gene elements, such as exon, intron, mrna, gene, splicing site, 
5'-flanking, 3'-flanking. Both cds and gene can be used in combination with 5'-flanking and 3'-flanking.

The position file can be psl or gff format now. The genome sequence file must
be fasta format. Note that this program is originally designed for CDS region, in other words,
not consider the UTR regions. You should alter on this, in order not to make mistake. 
For psl format, it only considers about blocks. For gff format, it only recognize mRNA and CDS feature.
The program will detect the file format automatically, or you can set by "--posformat" manually.

The default for "--type" option is "mrna", here it has the same meaning as cds, because we don't consider
UTR region.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2007-12-26

=head1 Usage
  % $0  [option] <pos_file> <seq_file>
  --posformat <str>   specify position file format
  --type <str>        specify element type: exon,intron,mrna,gene,splice,flank5,flank3, default=mrna
  --verbose           output verbose information to screen  
  --help              output help information to screen  

=head1 Exmple

 perl ./getGene.pl chr01.psl chr01.fa 
 perl ./getGene.pl chr01.psl chr01.fa -type mrna 

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Posformat,$Type,$Flank5,$Flank3);
my ($Verbose,$Help);
GetOptions(
	"type:s"=>\$Type,
	"flank5:i"=>\$Flank5,
	"flank3:i"=>\$Flank3,
	"posformat:s"=>\$Posformat,
	"verbose!"=>\$Verbose,
	"help!"=>\$Help
);
$Type ||= "mrna";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $pos_file = shift;
my $seq_file = shift;

my %gene;

read_gff($pos_file,\%gene) if($Posformat eq 'gff' || $pos_file =~ /.gff$/ || $pos_file =~ /\.gff\.gz$/);

#open(IN,$seq_file)||die("failed $seq_file\n");
if ($seq_file =~ /\.gz$/) {
	open IN, "gunzip -c $seq_file | ";
} else {
	open IN, $seq_file;
}

$/=">"; <IN>; $/="\n";	
while (<IN>) {
	my $output;
	my $chr=$1 if(/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\s//g;
	$/="\n";

	my $seq_len=length($seq);
	my $chr_pp=$gene{$chr};
	foreach  my $gene (sort keys %$chr_pp) {
		my $strand=$$chr_pp{$gene}{strand};
		next if(!exists $chr_pp->{$gene}{exon});
		my @exon = @{$chr_pp->{$gene}{exon}};

		if ($Type eq "mrna") {
			my $mrna;
			my ($left_leng, $right_leng);
			$left_leng = $exon[0][0]-1 if($left_leng > $exon[0][0]-1);
			$right_leng = $seq_len - $exon[-1][1] if($right_leng > $seq_len - $exon[-1][1]);

			$mrna .= substr($seq,$exon[0][0]-$left_leng-1,$left_leng) if($left_leng);
			for (my $i=0; $i<@exon; $i++) {
				$mrna .= substr($seq,$exon[$i][0]-1, $exon[$i][1] - $exon[$i][0] + 1);	

			}
			$mrna .= substr($seq,$exon[-1][1],$right_leng) if($right_leng);
			$mrna = Complement_Reverse($mrna) if($strand eq '-');
			Display_seq(\$mrna);
			my $mark = "$gene  [mRNA]";
			$output .= ">".$mark."\n".$mrna;
		}
	
		
	}
	print $output;
}

close(IN);




####################################################
################### Sub Routines ###################
####################################################

#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################


#############################################
sub Complement_Reverse{
	my $seq=shift;
	$seq=~tr/AGCTagct/TCGAtcga/;
	$seq=reverse($seq);
	return $seq;

}
#############################################



sub read_gff{
	my $file=shift;
	my $ref=shift;
	if ($file =~ /\.gz$/){
		open(IN, "zcat $file|") or die $!;
	}
	else{
		open (IN,$file) || die ("fail open $file\n");
	}
	while (<IN>) {
		s/^\s+//;
		s/\s+$//;
		my @t = split(/\s+/);
		my $tname = $t[0];
		my $qname;
		if ($t[2] eq 'mRNA' || $t[2] eq 'CDS' || $t[2] eq 'transcript') {
			$qname = $1 if($t[8] =~ /^GenePrediction\s+(\S+)/ || $t[8] =~ /^ID=([^;]+);*/ || $t[8] =~ /^Parent=([^;]+);*/);
			##print $qname."\n";
		}
		if ($t[2] eq 'match' || $t[2] eq 'HSP') {
			$qname = $1 if($t[8] =~ /Target\s+\"(\S+)\"/);
		}

		
		if ($t[2] eq 'mRNA' || $t[2] eq 'match' || $t[2] eq 'transcript') {
			$ref->{$tname}{$qname}{strand} = $t[6];
		}
		if ($t[2] eq 'CDS' || $t[2] eq 'HSP') {
			push @{$ref->{$tname}{$qname}{exon}}, [$t[3],$t[4]];
		}
	}
	close(IN);

	##print Dumper $ref;
	
	##change the exon order
	foreach my $chr (keys %$ref) {
		my $chr_p = $ref->{$chr};
		foreach my $gene (keys %$chr_p) {
			my $gene_p = $chr_p->{$gene};
			next if(!exists $gene_p->{exon});
			my @exon = sort {$a->[0] <=> $b->[0]} @{$gene_p->{exon}};
			$gene_p->{exon} = \@exon;
		}
	}
}
