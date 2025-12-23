#!/usr/bin/perl

=head1 Description

        This pipeline is to use kgf to fill gaps. kgf is a new gap filling program, based on PE reads and gap edge related reads. All the reads assembled by overlap method. This pipeline contains collecting gap related reads, reads filter, gap filling gap analysis, filling gaps  and filling result analysis.
	To get gap related reads, there are two choise. You can use krs to map sequencing reads to raw scaffold sequences, or use SOAPdenovo reads mapping result. This pipeline is used to get reads from SOAPdenovo map results.

=head1 Contact and Version
        Contact: Assemble Development Team(ADT).
        Version: 2.0
        Date: 2011.10.11

=head1 Usage
  
  perl SRKGF.pl [options]
  --dir		<str>	set SOAPdenovo working directory, reading file: *.scaf, *.scafSeq, *.shortreadInGap(.gz), *.PEreadOnContig(.gz), *.RlongReadInGap(.gz) , *contig for step 1
  --kmer	<int>	set kmer size used in SOAPdenovo. default 31 .  for step 1
  --prefix	<str>	set prefix used by SOAPdenovo for step 1
  --scaf	<str>	set scaffold sequence. for step 1
  --contig	<str>	set contig sequence. for step 2
  --gapread	<str>	set gap reads. for step 2
  --step	<int>	set steps. 1: get gap read from grape assemble result.2: contig gap filling 3: fill result analysis. default=123.
  --thread	<int>	set thread number for kgf gap filling, default =8. for step 2.
  --cpu		<int>	set scaffolds cut number for kgf gapfilling step, default=1. for step 1,2.
  --outdir	<str>	set the output directory. for all.
  --cvg		<int>	set the average coverage depth of contig . for step 1.
  --queue	<str>	set the queue for qsub jobs . for step 2.
  --vf		<str>	set vf for kgf.sh. for step 2.
  --shortcontig	<int>	set contig length cutoff, when contigs less than it, the contigs will be masked into gap. default 100 - (kmersize). for step 1.
  --falsecontig	<int>	when contigs less than it and related gap length is One, the shorter contig will be masked into gap. default 150.for step 1. 
  --nq 	      <local/qsub> switch option , no jobs will be qsub when set this option . default qsub . 	
  --verbose   output verbose information to screen  
  --help      output help information to screen

=head1 example
perl path/SRKGF.pl --dir path/SOAPdenovo_output_dir  --outdir path/SRKGF_output_dir --prefix SOAPdenovo_file_prefix --kmer SOAPdenovo_Kmer_set --cvg contig_coverage --cpu directory_num_for_kgf --thread kgf_thread  --step 123 .  --queue  compute_node.q  >logfile

=head1 
perl ./SRKGF.pl --dir /path/02.assemble/  --outdir /path/gap_fill/SRkgf/  --prefix test --kmer 31 --cvg 25 --cpu 8 --thread 8 --step 123 . --queue all.q  >fill.log 

=head1 
more example see Pipe_dir/example.sh . 

=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);

my ($path,$path2,$kgf2) = ("$Bin/stat2","$Bin/SR2","$Bin/kgf1.19");

my ($scafSeq,$lib_file,$contig_file,$gapread_file,$dir,$prefix,$IsGz,$noqsub,$Queue);
my ($cpu,$step,$iter, $Help, $Outdir, $Verbose, $thread, $startid,$scafnum, $kvf, $len, $vf, $kmer,$cvg,$contigcutoff1,$contigcutoff2,$longread);
my $minover;
GetOptions(
	"scaf:s"=>\$scafSeq,
#	"lib:s"=>\$lib_file,
	"dir:s"=>\$dir,
	"prefix:s"=>\$prefix,
	"contig:s"=>\$contig_file,
	"gapread:s"=>\$gapread_file,
	"kmer:i"=>\$kmer,
	"cvg:i"=>\$cvg,
	"cpu:i"=>\$cpu,
	"vf:s"=>\$kvf,
	"step:i"=>\$step,
#	"iter:i"=>\$iter,
	"thread:i"=>\$thread,
#	"lrs:i"=>\$longread,
	"shortcontig:i"=>\$contigcutoff1,
	"falsecontig:i"=>\$contigcutoff2,
	"outdir:s"=>\$Outdir,
	"queue:s"=>\$Queue,
	"nq!"=>\$noqsub,
	"verbose"=>\$Verbose,
	"help"=>\$Help,
	"minover"=>\$minover
);

$Outdir ||= ".";
$cpu ||= 4;
$step ||=123;
$iter ||=0;
$thread ||=4;
$startid ||=0;
$kvf ||="0G";
$scafnum ||=0;
#$longread ||=0;
$kmer ||=31;
$cvg ||=25;
$contigcutoff2 ||= 150;
$contigcutoff1 ||= 100 - $kmer;
$Queue ||= "";
#$IsGz ||= 0; #2011-4-18 !!
$minover||=30;
my ($scaf,$reads,$PEreads,$Rlongread,$totalread,$ctgfile);

#die `pod2text $0` if(!$prefix);
die `pod2text $0` if ($Help);

#2011-5-26 
sub abs_path
{
	chomp(my $tem_dir = `pwd`);
	foreach(@_){
		$_ || next ;
		/^\// && next ;
		$_ = "$tem_dir/$_";
		s/\/+$//;
	}
}
abs_path($Outdir,$dir,$contig_file,$gapread_file,$scafSeq);
#$Outdir =~ s/\/$//; #delete the last /
mkdir $Outdir unless (-d $Outdir);

my $size=0;
my $Rlongflag = 0 ;
#my $pwd = `pwd`;
#if($Outdir == "\.")
#{
#	$Outdir = $pwd;
#}

print "Starting kgf gap filling , with step: $step \n";


if($step =~ /1/)
{
	die `pod2text $0` unless($dir);
	die `pod2text $0` if(!$prefix);
	print "Step one: get reads related with gaps from the result of SOAPdenovo!\n";
	$dir =~ s/\/$//;
	$scaf = $dir."/$prefix.scaf";
	$scafSeq = $dir."/$prefix.scafSeq";
	
	#ensure all necessary files exist.
	unless (-e $scaf && -e $scafSeq)
	{
		print "check file :\n$scaf\nand file :\n$scafSeq\n";
	}

	$reads = $dir."/$prefix.shortreadInGap";
	unless(-e $reads){
		$reads = $dir."/$prefix.shortreadInGap.gz";
		unless(-e $reads)
		{
			print "check file :\n$reads\n";
		}
	}

	$PEreads = $dir."/$prefix.PEreadOnContig";
	unless (-e $PEreads){
		$PEreads = $dir."/$prefix.PEreadOnContig.gz";
		unless (-e $PEreads)
		{
			print "check file :\n$PEreads\n"
		}
	}

	$Rlongread = $dir."/$prefix.RlongReadInGap";
	unless (-e $Rlongread){
		$Rlongread = $dir."/$prefix.RlongReadInGap.gz";
		unless (-e $Rlongread)
		{
			print "check file :\n$Rlongread\n";
			$Rlongflag = 1 ;
		}
	}

	$ctgfile = $dir."/$prefix.contig";

	my $linkscaf = "$Outdir/$prefix.scafSeq";
	if(-e $linkscaf)
	{
		$scafSeq = $linkscaf;
	}else{
		`ln -s $scafSeq $Outdir/$prefix.scafSeq`;
		$scafSeq = $linkscaf;
	}

	if($contigcutoff1 < 70)
	{
		$contigcutoff1 = 70;
	}
	my $n = int($contigcutoff1 / $kmer);
	if($n > 3) { $n = 3 ;}
	print "getreadgap use the offset $n * $kmer\n" ;
	$kmer  = $n * $kmer ;
	`$path2/ggi -s $scafSeq -f $scaf -c $contigcutoff1 -t $contigcutoff2 -o $Outdir/$prefix.scafSeq -p $thread >$Outdir/gapinfo.log 2>$Outdir/ggi.error`;
	
	#ggr has a bug which the *.PEreadOncontig.gz must be generated by SOAPdenovo2 , otherwise ggr will be broken down .
	open WS , ">$Outdir/ggr.sh" ;
	#print WS "$path2/ggr -g $Outdir/$prefix.scafSeq.gapInfo -s $reads  -p $PEreads -K $kmer -d $cvg -o $Outdir -m $Outdir/$prefix.scafSeq.maskContig -c $ctgfile >$Outdir/getgapread.log 2>$Outdir/getgapread.error\n" ;
	if($Rlongflag){
	print WS "$path2/ggr -g $Outdir/$prefix.scafSeq.gapInfo -s $reads -p $PEreads -K $kmer -d $cvg -o $Outdir -m $Outdir/$prefix.scafSeq.maskContig -c $ctgfile >$Outdir/getgapread.log 2>$Outdir/getgapread.error\n" ;
    	}else{
	print WS "$path2/ggr -g $Outdir/$prefix.scafSeq.gapInfo -s $reads -r $Rlongread -p $PEreads -K $kmer -d $cvg -o $Outdir -m $Outdir/$prefix.scafSeq.maskContig -c $ctgfile >$Outdir/getgapread.log 2>$Outdir/getgapread.error\n" ;
	}
	my $gapInfoFile = "$Outdir/$prefix.scafSeq.gapInfo" ;
	my $maskCtgFile = "$Outdir/$prefix.scafSeq.maskContig" ;
	my $gvf = (-s $gapInfoFile) + (-s $maskCtgFile) + (-s $ctgfile);
	$gvf = $gvf/1000000 ;
	if($gvf > 1000) #1G
	{
		$gvf = $gvf/1000;
		$gvf = $gvf."G";
	}else {
		$gvf = $gvf."M" ;
	}

	print "ggr use memory : $gvf\n" ;
	if($noqsub){
		`sh $Outdir/ggr.sh` ;
	}else{
	    	if($Queue eq ""){
			`nohup $path/qsub-sge.pl --resource vf=$gvf --maxjob 1 --jobprefix gr --convert no $Outdir/ggr.sh`;
		}else{
		`nohup $path/qsub-sge.pl --resource vf=$gvf --maxjob 1 --queue $Queue --jobprefix gr --convert no $Outdir/ggr.sh` ;
	    	}
	}
	#`$path2/ggr -g $Outdir/$prefix.scafSeq.gapInfo -s $reads -r $Rlongread -p $PEreads -K $kmer -d $cvg -o $Outdir/gapread.fa -m $Outdir/$prefix.scafSeq.maskContig -c $ctgfile >$Outdir/getgapread.log 2>$Outdir/getgapread.error`;
	$contig_file = "$scafSeq.SCAF.contig";
	if(-e "$Outdir/shortread.fa.gz" && -e "$Outdir/PEread.fa.gz" && -e "$Outdir/longread.fa.gz")
	{
		`cat $Outdir/shortread.fa.gz $Outdir/PEread.fa.gz $Outdir/longread.fa.gz >$Outdir/gapread.fa.gz` ;
		`rm $Outdir/shortread.fa.gz $Outdir/PEread.fa.gz $Outdir/longread.fa.gz`;
	}elsif(-e "$Outdir/shortread.fa.gz" && -e "$Outdir/PEread.fa.gz"){
		`cat $Outdir/shortread.fa.gz $Outdir/PEread.fa.gz >$Outdir/gapread.fa.gz`;
		`rm $Outdir/shortread.fa.gz $Outdir/PEread.fa.gz` ;
	}
	#$gapread_file = "$Outdir/gapread.fa";
	#`gzip $gapread_file`;
	$gapread_file = "$Outdir/gapread.fa.gz";
	unless(-e $gapread_file){
		print "gzip error , there has no gapread.fa.gz file\n";
	}
	`perl $path/gz_readstat2.pl $gapread_file >$Outdir/gapread.fa.gz.stat`;
	
	print "Get gap read finished!\n";	
}


if($step=~/2/ && $cpu > 1)
{
	print "step two: begin to Cut the files and kgf gap filling!";
	if(!$contig_file)
	{
		if($scafSeq){  #2011-4-19
			my @names=split (/\//,$scafSeq);
			my $name = $names[@names -1];
			$scafSeq = "$Outdir/$name";
		}
		if(!$scafSeq){
			print "please set scaffold sequence!\n";
			exit(1);
		}
		my $SCAFctg = "$scafSeq.SCAF.contig";
		unless (-e $SCAFctg)
		{
			print "please make sure the progrom can find the file $scafSeq.SCAF.contig in this work directory!\n";
			exit(1);
		}
		$contig_file = "$scafSeq.SCAF.contig";
	}
	
	$gapread_file = "$Outdir/gapread.fa.gz" if(!$gapread_file);

	if(-e $gapread_file && -e $contig_file)
	{
		print "Load gap file: $gapread_file.\nLoad contig file: $contig_file.\n";
	}else{
		print "Please input or check $gapread_file or $contig_file!\n";
	}
	
	`perl $path/gz_Cut2.pl $gapread_file $contig_file $cpu $Outdir`;
	`perl $path/gz_Creatkgf.pl $kgf2 $Outdir $cpu $thread $minover >$Outdir/kgf.sh`;

	if($kvf eq '0G')
	{
	    	#$kvf=$thread*2;# + int($len/1000000000/$cpu);
		my $gapread_size = -s $gapread_file ;
		#$gapread_size *= 4 if($gapread_file =~ /(\w+)\.gz/);
		my $contig_size = -s $contig_file ;
		print "Calculate kgf necessary  memory : size(gapread.fa)/(--cpu) + size(*.SCAF.contig)*2/(--cpu)\n\n" ;
		$kvf = ($gapread_size/$cpu/1000000000) + ($contig_size/$cpu/1000000000)*2 ;
	    	$kvf=$kvf."G";
	}
	print "begin to fill gaps with kvf: $kvf\n";
	if($noqsub){
		`sh $Outdir/kgf.sh`;
	}else{
	    	if($Queue eq ""){
			`nohup $path/qsub-sge.pl --resource vf=$kvf --maxjob $cpu --jobprefix fg --convert no $Outdir/kgf.sh`;
		}elsif($Queue =~ /(w+)\.q$/){
			`nohup $path/qsub-sge.pl --resource vf=$kvf --queue $Queue --maxjob $cpu --jobprefix fg --convert no $Outdir/kgf.sh`;
	    	}else{
			`nohup $path/qsub-sge.pl --resource vf=$kvf --maxjob $cpu --jobprefix fg --convert no $Outdir/kgf.sh` ;
		}
	}
	`cat $Outdir/F*/FilledScaf/seq.thread* $scafSeq.CONTIG >$scafSeq.fill`;
	`cat $Outdir/F*/Log/log.thread* >$Outdir/fill.Log`;
	`cat $Outdir/F*/Snp/snp.thread* >$Outdir/fill.snp`;

	my $seq = "$Outdir/gapSeq.fa";
	`cat $Outdir/F*/gapSeq.fa > $seq`;

	`perl $path/get_scaftig.pl $scafSeq.fill >$scafSeq.fill.scaftig`;
	`perl $path/seq_n50 $scafSeq.fill.scaftig >$Outdir/fill.fasta.0.scaftig.N50`;
	my $n50 = `grep 'N50' $Outdir/fill.fasta.0.scaftig.N50 |awk '{print $2}'`;

	my $fill = `perl $path/gapfillratio2.pl $Outdir $cpu `;
	if($fill=~ /Error/)
	{
		print "\n\n$fill\n\n";
		exit(1);
	}								              
	print $fill;
	print "Gap filling finished! \n";

}

if($step =~ /3/)
{
	print "Step Three: analysis the fill result!\n";
	#2011-4-19
	my $depth = "$Outdir/gapread.fa.gz.gapread.depth";
	my $seq = "$Outdir/gapSeq.fa";
	my $Log = "$Outdir/fill.Log";
	my @info=glob("$Outdir/*.scafSeq.gapInfo");
	my $gapinfo = $info[0];
	unless (-e $depth && -e $seq && -e $Log && -e $gapinfo)
	{
		print "please check file: $depth, $seq , $Log and $gapinfo\n";
		print "make sure the file $depth , $seq , $Log and *.scafSeq.gapInfo under the dir to output!\n";
		exit(1);
	}
	`grep 'TRGAP' $Log >$Outdir/TR.lst`;
	`grep 'ERGAP' $Log >$Outdir/ER.lst`;
	`perl $path/fullfill.pl $seq >$Outdir/fullfill.lst`;
	`perl $path/get_unfull_info.pl $Outdir/fullfill.lst $gapinfo >$Outdir/unfill.gapInfo 2>$Outdir/unfill.stat`;

	`perl $path/blank.pl $depth $seq >$Outdir/blank.lst`;
	`perl $path/covercheck.pl $Outdir/blank.lst >$Outdir/cover.lst`;
	#`more $Outdir/cover.lst |awk '$3 < 1 {print $_}' >$Outdir/uncover.lst`;
	`perl $path/gether.pl $Outdir/cover.lst $Outdir/TR.lst $depth $Outdir/ER.lst $Outdir >$Outdir/canfill.stat`;

	print "Step Three finished! Please read file: unfill.stat and canfill.stat \n";
	#`rm $Outdir/gapread.fa`;
}


print "All Pipeline finished!";
