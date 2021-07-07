#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin $Script);
die "Usage: <config> <outdir> <priority[1:intact ORF; 2:longest pep length]\n" unless @ARGV == 3;
my $config = shift;
my $outdir = shift;
my $priority = shift;
mkdir $outdir unless -e $outdir;

my $qsub = 1;

## scritp needed
my $gtf2gff = "$Bin/gtf2gff.pl";
my $renameGenome = "$Bin/renameGenome.pl";
my $check_orf_for_gff = "$Bin/check_orf_for_gff.pl";
my $bestProteinFromEnsembl = "$Bin/bestProteinFromEnsembl.pl";
my $select_gff = "$Bin/select_gff.pl";
my $getGene = "$Bin/getGene.pl";
my $cds2aa = "$Bin/cds2aa.pl";
my $select_orf = "$Bin/select_orf.pl";
my $check_orf_for_cds = "$Bin/check_orf_for_cds.pl";

foreach my $p (\$gtf2gff, \$renameGenome, \$check_orf_for_gff, \$bestProteinFromEnsembl, \$select_gff, \$getGene, \$cds2aa, \$select_orf, \$check_orf_for_cds) {
	die "$$p does not exist!\n" unless -e $$p;
}

my %speciesFiles;
&readConfig($config, \%speciesFiles);

my %cmd;
foreach my $sp (keys %speciesFiles) {
	push @{$cmd{$sp}}, "perl $renameGenome $speciesFiles{$sp}{fa} >$outdir/$sp.fa";
	push @{$cmd{$sp}}, "perl $gtf2gff $speciesFiles{$sp}{gtf} >$outdir/$sp.all.gff";
	push @{$cmd{$sp}}, "perl $check_orf_for_gff $outdir/$sp.all.gff $outdir/$sp.fa >$outdir/$sp.all.gff.orf";
	push @{$cmd{$sp}}, "perl $bestProteinFromEnsembl $speciesFiles{$sp}{pep} $outdir/$sp.all.gff.orf $priority >$outdir/$sp.pep.best";
	push @{$cmd{$sp}}, "perl $select_gff $outdir/$sp.pep.best $outdir/$sp.all.gff >$outdir/$sp.best.gff";
	push @{$cmd{$sp}}, "perl $getGene $outdir/$sp.best.gff $outdir/$sp.fa >$outdir/$sp.best.cds";
	#push @{$cmd{$sp}}, "perl $cds2aa --check $outdir/$sp.best.cds | awk '\$4==0' >$outdir/$sp.best.cds.preStop";
	push @{$cmd{$sp}}, "perl $check_orf_for_cds $outdir/$sp.best.cds | awk '\$4>0' >$outdir/$sp.best.cds.preStop";
	push @{$cmd{$sp}}, "perl $cds2aa $outdir/$sp.best.cds >$outdir/$sp.best.pep";
	push @{$cmd{$sp}}, "perl $select_orf $outdir/$sp.best.pep $outdir/$sp.all.gff.orf >$outdir/$sp.best.gff.orf";
	push @{$cmd{$sp}}, "awk '\$8==1 && \$9==1' $outdir/$sp.best.gff.orf >$outdir/$sp.best.gff.orf.intact";
	push @{$cmd{$sp}}, "rm $outdir/$sp.all.gff $outdir/$sp.all.gff.orf $outdir/$sp.pep.best";
}

my @shellFiles;
foreach my $sp (keys %cmd) {
	my $sh_file = "$outdir/$sp.sh";
	open SH, ">$sh_file";
	print SH "date\n";
	foreach my $line (@{$cmd{$sp}}) {
		print SH "$line\n";
	}
	print SH "date\n";
	close SH;
	push @shellFiles, $sh_file;
}
foreach my $sh_file (@shellFiles) {
	if ($qsub) {
		#system "qsub -S /bin/sh -cwd -l vf=0.01G -q st_pc.q -P pc_pa_un $sh_file";
		system "qsub -S /bin/sh -cwd -l vf=0.01G -q ngb.q -P ngb_un $sh_file";
	} else {
		system "sh $sh_file";
	}
}

################################
sub readConfig {
	my ($in_file, $ref) = @_;
	open IN, $in_file;
	$/ = "##";
	<IN>;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//g;
		my @lines = (split /\s+/);
		die unless @lines == 4;
		my $sp = $lines[0];
		$ref->{$sp}{fa} = $lines[1];
		$ref->{$sp}{pep} = $lines[2];
		$ref->{$sp}{gtf} = $lines[3];
	}
	close IN;
}
