#!/usr/bin/perl -w
use strict;
die "Usage: <in_file(gff|psl)> <file type> <personal region>\n" unless @ARGV >= 2;
my $file_type = $ARGV[1];

## get personal information.
my %personal = ();
&get_personalInfo($ARGV[2], \%personal) if $ARGV[2];

my (%cdsPos, %utrPos);
&get_exonPos($ARGV[0], \%cdsPos, \%utrPos, $file_type);
my %pos;
foreach my $gene (keys %cdsPos) {
	foreach my $p (@{$cdsPos{$gene}}) {
		my ($chr, $bg, $ed, $strand) = @$p; 
		push @{$pos{$gene}}, ($bg, $ed);
	}
}
foreach my $gene (keys %utrPos) {
	foreach my $p (@{$utrPos{$gene}}) {
		my ($chr, $bg, $ed, $strand) = @$p; 
		push @{$pos{$gene}}, ($bg, $ed);
	}
}

my $widest = 0;
my $geneCount = 0;
my ($region_bg, $region_ed);
foreach my $gene (keys %pos) {
	$geneCount ++;
	@{$pos{$gene}} = sort {$a <=> $b} @{$pos{$gene}};
	my $len = $pos{$gene}[-1] - $pos{$gene}[0] + 1;
#$widest = $len if $widest < $len;
	if ($region_bg) {
		$region_bg = $pos{$gene}[0] if $pos{$gene}[0] < $region_bg;
	} else {
		$region_bg = $pos{$gene}[0];
	}
	if ($region_ed) {
		$region_ed = $pos{$gene}[-1] if $pos{$gene}[-1] > $region_ed;
	} else {
		$region_ed = $pos{$gene}[-1];
	}
	$widest = $region_ed - $region_bg + 1;
}
my $rate = $widest/700;
my $height = $geneCount * 80;


print  '<?xml version="1.0" standalone="no"?>', "\n";
print  '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">', "\n";
print  "<svg width=\"900\" height=\"$height\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">", "\n";
print "<rect x=\"0\" y=\"0\" width=\"900\" height=\"$height\" fill=\"white\"/>\n";

my ($leftBar, $topBar) = (50, 50);
my ($line_y1) = ($topBar);
foreach my $gene (keys %cdsPos) {
	my ($chr_bg, $chr_ed) = ($pos{$gene}[0], $pos{$gene}[-1]);
	my $region_len = $pos{$gene}[-1] - $pos{$gene}[0] + 1;
	my $line_len = $region_len/$rate;
	my $line_x1 = ($chr_bg - $region_bg)/$rate +  $leftBar - 20;
	my $line_x2 = $line_x1 + $line_len + 20;
	my $line_y2 = $line_y1;
	printf '<text x="%d" y="%d" font-family="Arial" font-size="12" fill="black">%s</text>' . "\n", $line_x1, $line_y1-25, $gene; 
	print "<line x1=\"$line_x1\" y1=\"$line_y1\" x2=\"$line_x2\" y2=\"$line_y2\" style=\"stroke:rgb(99,99,99);stroke-width:1\"/>\n";

	my $count = 0;
	foreach my $p (@{$cdsPos{$gene}}) {
		$count ++;
		my ($chr, $bg, $ed, $strand) = @$p;
		my $block_color = ($chr =~ /yellow/i) ? "yellow" : "red";

		my $cds_len = $ed - $bg + 1;
		my $rect_len = $cds_len/$rate;
#print "test:$cds_len\t$rate\n";
		my $rect_x  = $leftBar + ($bg - $region_bg + 1)/$rate;
		my $rect_y = $line_y1 - 5;
		my $text_x = $rect_x - 5;
#my $text_y = $count % 2 ? $rect_y - 10 : $rect_y + 30;
		my $text_y = $rect_y + 30;
		printf '<text x="%d" y="%d" font-family="Arial" font-size="8" fill="%s">%s</text>' . "\n", $text_x, $text_y, ,"black", $chr if $count == 1;

		print "<rect x=\"$rect_x\" y=\"$rect_y\" width=\"$rect_len\" height=\"10\" style=\"fill:$block_color;stroke:$block_color;stroke-width:0.5; fill-opacity:1;stroke-opacity:1\"/>\n";

		my ($tri_y1, $tri_y2, $tri_y3) = ($line_y1, $line_y1+10, $line_y1-10);
		my ($tri_x1, $tri_x2, $tri_x3);
		if ($strand eq "+") {
			$tri_x1 = $rect_x + $rect_len + 15;
			$tri_x2 = $rect_x + $rect_len;
			$tri_x3 = $rect_x + $rect_len;
		} elsif ($strand eq "-") {
			$tri_x1 = $rect_x - 15;
			$tri_x2 = $rect_x;
			$tri_x3 = $rect_x;
		}
		if ($strand eq "+") {
			print "<polygon points=\"$tri_x1,$tri_y1 $tri_x2,$tri_y2 $tri_x3,$tri_y3\" style=\"fill:$block_color; stroke:$block_color;stroke-width:0.5\"/>\n" if $count == @{$cdsPos{$gene}};
		} else {
			print "<polygon points=\"$tri_x1,$tri_y1 $tri_x2,$tri_y2 $tri_x3,$tri_y3\" style=\"fill:$block_color; stroke:$block_color;stroke-width:0.5\"/>\n" if $count == 1;
		}
	}
=pod
	$count = 0;
	foreach my $p (@{$utrPos{$gene}}) {
		my ($chr, $bg, $ed, $strand) = @$p;
		my $block_color = ($chr =~ /yellow/i) ? "yellow" : "white";

		my $cds_len = $ed - $bg + 1;
		my $rect_len = $cds_len/$rate;
#print "test:$cds_len\t$rate\n";
		my $rect_x  = $leftBar + ($bg - $region_bg + 1)/$rate;
		my $rect_y = $line_y1 - 5;
		my $text_x = $rect_x - 5;
#my $text_y = $count % 2 ? $rect_y - 10 : $rect_y + 30;
		my $text_y = $rect_y + 30;
		printf '<text x="%d" y="%d" font-family="Arial" font-size="8" fill="%s">%s</text>' . "\n", $text_x, $text_y, ,"black", $chr if $count == 1;

		print "<rect x=\"$rect_x\" y=\"$rect_y\" width=\"$rect_len\" height=\"10\" style=\"fill:$block_color;stroke:black;stroke-width:0.5; fill-opacity:1;stroke-opacity:1\"/>\n";

		my ($tri_y1, $tri_y2, $tri_y3) = ($line_y1, $line_y1+10, $line_y1-10);
		my ($tri_x1, $tri_x2, $tri_x3);
		if ($strand eq "+") {
			$tri_x1 = $rect_x + $rect_len + 15;
			$tri_x2 = $rect_x + $rect_len;
			$tri_x3 = $rect_x + $rect_len;
		} elsif ($strand eq "-") {
			$tri_x1 = $rect_x - 15;
			$tri_x2 = $rect_x;
			$tri_x3 = $rect_x;
		}
		if ($strand eq "+") {
			print "<polygon points=\"$tri_x1,$tri_y1 $tri_x2,$tri_y2 $tri_x3,$tri_y3\" style=\"fill:$block_color; stroke:$block_color;stroke-width:0.5\"/>\n" if $count == @{$cdsPos{$gene}};
		} else {
			print "<polygon points=\"$tri_x1,$tri_y1 $tri_x2,$tri_y2 $tri_x3,$tri_y3\" style=\"fill:$block_color; stroke:$block_color;stroke-width:0.5\"/>\n" if $count == 1;
		}
	}
=cut


	## draw personal regions.
	foreach my $p (@{$personal{$gene}}) {
		my ($chr, $bg, $ed, $color) = @$p;
		my $rect_len = ($ed - $bg + 1)/$rate;
		my $rect_x  = $leftBar + ($bg - $region_bg + 1)/$rate;
		my $rect_y = $line_y1 - 150;
		my $block_color = $color;
		my $stroke_color = $color;
		print "<rect x=\"$rect_x\" y=\"$rect_y\" width=\"$rect_len\" height=\"300\" style=\"fill:$block_color;stroke:$stroke_color;stroke-width:1; fill-opacity:1;stroke-opacity:1\"/>\n";
	}


	$line_y1 += 80;
}


print  '</svg>', "\n";


## subroutine
#######################
sub get_chrGene {
	my $in_file = shift;
	my $ref1 = shift;
	my $ref2 = shift;
	my $ref3 = shift;
	open IN, $in_file;
	while (<IN>) {
		chomp;  
		my @info = split /\t/;
		next unless $info[2] eq 'mRNA';
		$info[-1] =~ /ID=(.+?);/;
		my $gene = $1; 
		push @{$ref1->{$info[0]}}, [$gene, $info[3], $info[4], $info[6]];
	}
	close IN;
	foreach my $chr (keys %$ref1) {
		my %p_to_bg;
		foreach my $p (@{$ref1->{$chr}}) {
			my $bg = $p->[1];
			$p_to_bg{$p} = $bg;
		}
		@{$ref1->{$chr}} = sort {$p_to_bg{$a} <=> $p_to_bg{$b}} @{$ref1->{$chr}};
		for (my $i = 0; $i < @{$ref1->{$chr}}; $i ++) {
			my $p = $ref1->{$chr}->[$i];
			my ($gene, $bg, $ed, $strand) = @$p;
			my $index = $i + 1;
			push @$p, $index;
			$ref2->{$gene} = [$chr, $bg, $ed, $strand, $index];
			$ref3->{$chr}{$index} = $gene;
		}
	}
}

## get cds pos of each gene
sub get_exonPos {
	my ($in_file, $cds_ref, $utr_ref, $file_type) = @_;
	open IN, $in_file;
	if ($file_type eq "gff") {
		while (<IN>) {
			my @info = split /\s+/;
			if ($info[2] eq "CDS") {
				#next unless $info[2] eq "CDS";
				die unless $info[8] =~ /Parent=(\S+?);/ || $_ =~ /gene_id\s+"(\S+?)";/;
				my $id = $1;
				push @{$cds_ref->{$id}}, [$info[0], $info[3], $info[4], $info[6]];
			} elsif ($info[2] =~ /UTR/) {
				die unless $info[8] =~ /Parent=(\S+?);/;
				my $id = $1;
				push @{$utr_ref->{$id}}, [$info[0], $info[3], $info[4], $info[6]];
			}
#push @$ref2, $id;
		}
	} elsif ($file_type eq "psl") {
		while (<IN>) {
			my @info = split /\s+/;
			my $id = $info[9];
			my $chr = $info[13];
			my $strand = $info[8];
			my @block_len = split /,/, $info[-3];
			my @query_pos = split /,/, $info[-2];
			my @refer_pos = split /,/, $info[-1];
			for (my $i = 0; $i < @block_len; $i ++) {
				my $len = $block_len[$i];
				my $q_bg = $query_pos[$i] + 1;
				my $q_ed = $q_bg + $len - 1;
				my $r_bg = $refer_pos[$i] + 1;
				my $r_ed = $r_bg + $len - 1;
				push @{$cds_ref->{$id}}, [$chr, $r_bg, $r_ed, $strand];
			}
		}
	}
	close IN;
	foreach my $id (keys %$cds_ref) {
		@{$cds_ref->{$id}} = sort {$a->[1] <=> $b->[1]} @{$cds_ref->{$id}};
	}
}

####################################################################
## store the personal information
sub get_personalInfo {
	my ($in_file, $ref) = @_;
	open IN, $in_file;
	while (<IN>) {
		next if /^#/;
		my @info = split /\s+/;
		push @{$ref->{$info[0]}}, [$info[1], $info[2], $info[3], $info[4]]; ## chr, bg, ed, color
	}
	close IN;
}
####################################################################

