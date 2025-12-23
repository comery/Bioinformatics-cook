#!/usr/bin/perl -w

=head1 Name
 GATK.pl   high-accuracy call SNP and indel with GATK walker.

=head1 Author

 Written by BGI-tech in april 2013
 
=head1 Flowchart
 
 
 CleanReads-----BWA----mark duplicate-----Local realignment----BQSR----Reduce Reads----UnifiedGenotyper call snp----VQSR----Analysis-ready Variants
 
=head1 Usage
 options:(write in main.sh )

  -genome *              genome file,must be end of .fa or .fasta
  -outdir *              output directory 
  -lib    *              sample info file 
  -queue                 specify queue
  -project               project ID 
  -help                  output help information to screen
  
=head1  Example of the main.sh:

  perl  GATK.pl -genome hg19.fasta -lib input.lib  -queue bc.q -project paeatest -outdir result

=cut

use strict;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
use FindBin '$Bin';
use lib $Bin;

my ($ProjectName,$LIB,$OUTDIR,$GENOME,$v,$k,$queue,$project,$Help);

GetOptions(
  "genome:s"=>\$GENOME,
	"lib:s"=>\$LIB,
	"outdir:s"=>\$OUTDIR,
	"queue:s" =>\$queue,
	"project:s" => \$project,
	"help|?" => \$Help
);
$queue = "bc.q" if (!defined $queue);  
$project = "paeatest"  if (!defined $project);
die `pod2text $0` if($Help||!defined$GENOME||!defined$LIB||!defined$OUTDIR);

my $dos2unix_dir = "$OUTDIR/data_dos2unix";
$LIB = abs_path($LIB);
$OUTDIR = abs_path($OUTDIR);

#*******************************************************************************************
# 初始化软件路径  initialize all dir used in pipeline
#*******************************************************************************************
my ($GATK,$BWA,$samtoolsdir,$samtools,$picard,$python,$BinDir);
if (-e '/ifs2') {
	$BinDir="/ifs2/BC_GAG/Group/zhuangzhenhua/bin/GATK"; #
	$GATK="java -jar /ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/GenomeAnalysisTK-2.4-9/GenomeAnalysisTK.jar"; #/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_CSAP/DNA_CSAP_5.2.4/bin/Tool/GenomeAnalysisTK.jar;
	$BWA="/ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/bwa-0.7.0/bwa"; #i/ifs2/PC_PA_UN_2/USER/zhouchengran/bin/bwa-0.7.10/bwa;
	$samtoolsdir="/ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/samtools-0.1.18"; #/ifs2/PC_PA_UN_2/USER/zhouchengran/bin/samtools-0.1.19;
	$samtools="/ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/samtools-0.1.18/samtools"; #/ifs2/PC_PA_UN_2/USER/zhouchengran/bin/samtools-0.1.19/samtools;
	$picard="/ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/picard/picard-tools-1.86"; 
	$python="/ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/Python-2.7/bin/python2.7";
}

#转成绝对路径 change to absolute directory
$GENOME = abs_path(dirname($GENOME)).'/'.basename($GENOME);

#*******************************************************************************************
# 读lib配置文件，获取样品信息
#*******************************************************************************************
my %data_info;
my ($KeyName,$FastqA,$FastqB,%KEYNAMES);
open LIB,"<$LIB" || die "can`t open the lib file:$!";
while (<LIB>)
{
	chomp($_);
	next if (/^\#/);
	next if (/^\s*$/);
my 	@info=split /\s+/,$_;
		$KeyName=$info[0];
		$FastqA = $info[1];
		die "can't find file: $FastqA" if(!-e $FastqA);
		$FastqA = abs_path($FastqA);
		$FastqB = $info[2] ;
		die "can't find file: $FastqB" if(!-e $FastqB);
		$FastqB = abs_path($FastqB);
		$KEYNAMES{$KeyName}->[0] = $FastqA;
		$KEYNAMES{$KeyName}->[1] = $FastqB;
		}
	close LIB;
	
my $INDEXDir="$OUTDIR/INDEX";
my $InFileDir="$OUTDIR/Inputfile";
my $BWADir="$InFileDir/bwa_reslt";
my $ChrbamDir="$BWADir/Chrbam";
my $PREDir="$InFileDir/pre_result";
my $PARTSHELLDir="$PREDir/shell";
my $CALLSNPDir="$OUTDIR/Callsnp";
my $CALLSNPLISTDir="$CALLSNPDir/list";
my $CALLSNPshellDir="$CALLSNPDir/shell";
my $CALLINDELDir="$OUTDIR/Callindel";
my $CALLINDELshellDir="$CALLINDELDir/shell";

MakeDir();


my @tissueArray=sort(keys%KEYNAMES); # 存放样品名称的数组 samples' names
my $TISSUENUM=keys(%KEYNAMES); # 样本数量 number of samples

my ($GENOMEID,@chr);
&showLog("cp genome file and generate script files,please waitting....");
open CLEAN, ">$OUTDIR/clean.sh" or die $!;
open INDEX , ">$INDEXDir/bulid_index.sh" || die  $!;
  if (defined $GENOME ) {
	system("cp $GENOME $INDEXDir");
	$GENOME="$INDEXDir/" . basename($GENOME);
	$GENOME=abs_path(dirname($GENOME)) . "/" . basename($GENOME);
	$GENOME=~s/(\S+)\.(\S+)$/$1/;
  $GENOMEID=$GENOME;
	$GENOME=$GENOME . ".fa";
#print $GENOME;
}
	open IN,$GENOME;
	while(<IN>){
		chomp;
		if(/\>(\S+)/){
			push @chr,$1;	#染色体ID
		}
	}
	close IN;
	print INDEX "$BWA index -a  bwtsw $GENOME  -p $GENOME\n";
	print INDEX "java -jar $picard/CreateSequenceDictionary.jar R=$GENOME  O=$GENOMEID.dict\n";
	print INDEX "$samtools  faidx $GENOME\n";

close INDEX;	
foreach my $k(@tissueArray)
{
	system("mkdir -m 755 -p $PREDir/$k") if (!-d "$PREDir/$k");
}

my ($MEANQUAL,$PARTSHELL);
open BWA, ">$BWADir/bwa.sh" or die $!;
open Chrbam, ">$ChrbamDir/Chrbam.sh" or die $!;
open QSUBPRE ,">$PARTSHELLDir/qsubpre.sh" or die $!;

 for my $k(sort keys %KEYNAMES) {
 	print BWA "$python $BinDir/bwaRunner.py PE  -v 0.04 -x 650 -l 35 -R 20 -t 4 -e 30 -i 15   $GENOME $KEYNAMES{$k}->[0] $KEYNAMES{$k}->[1] -d $BWADir  -s $k\n";
 	print Chrbam "$python $BinDir/sort_mergeRunner.py -t $samtools -p 4 -o $ChrbamDir  -s $k $BWADir/$k.sam.gz\n";
 	foreach my $j(@chr){
 		open  PARTSHELL ,">$PARTSHELLDir/$k.$j.sh" or die $!;
 		print PARTSHELL "#sort,mark duplicate,local Realigner\n";
 		print PARTSHELL "java  -jar   $picard/SortSam.jar  INPUT=$ChrbamDir/$k/$j/$k.$j.bam  OUTPUT=$PREDir/$k/$k.$j.sort.bam  SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT\n";
 		print PARTSHELL "java  -jar   $picard/AddOrReplaceReadGroups.jar  VALIDATION_STRINGENCY=LENIENT   I=$PREDir/$k/$k.$j.sort.bam O=$PREDir/$k/$k.$j.head.bam  ID=$k  LB=$k  PL=illumina SM=$k PU=$k\n";
 		print PARTSHELL "java  -jar   $picard/MarkDuplicates.jar MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000  INPUT=$PREDir/$k/$k.$j.head.bam OUTPUT=$PREDir/$k/$k.$j.dedup.bam  METRICS_FILE=$PREDir/$k/$k.$j.dedup.metrics  VALIDATION_STRINGENCY=SILENT\n";
 		print PARTSHELL "$samtools index $PREDir/$k/$k.$j.dedup.bam\n";
 		print PARTSHELL "$GATK -T RealignerTargetCreator -R $GENOME -I $PREDir/$k/$k.$j.dedup.bam -o  $PREDir/$k/$k.$j.realn.intervals\n";
 		print PARTSHELL "$GATK -T IndelRealigner -R $GENOME -I $PREDir/$k/$k.$j.dedup.bam -targetIntervals $PREDir/$k/$k.$j.realn.intervals -o $PREDir/$k/$k.$j.realign.bam\n";
 		print PARTSHELL "#generate confidence snp database file\n";
 		print PARTSHELL "$GATK  -T UnifiedGenotyper -stand_call_conf 30.0 -stand_emit_conf 0 -glm BOTH  -rf BadCigar  -R $GENOME -I $PREDir/$k/$k.$j.realign.bam -o $PREDir/$k/$k.$j.gatk.raw.vcf -nt 2 -nct 4\n";
 		print PARTSHELL "$samtools mpileup -ugf $GENOME  $PREDir/$k/$k.$j.realign.bam | $samtoolsdir/bcftools/bcftools view -Ncvg - > $PREDir/$k/$k.$j.samtools.raw.vcf\n";
 		print PARTSHELL "$GATK -T SelectVariants -R $GENOME --variant $PREDir/$k/$k.$j.gatk.raw.vcf --concordance $PREDir/$k/$k.$j.samtools.raw.vcf -o $PREDir/$k/$k.$j.concordance.raw.vcf\n";
 		print PARTSHELL "MEANQUAL=`awk '{ if (\$1 !~ /#/) { total += \$6; count++ } } END { print total/count }' $PREDir/$k/$k.$j.concordance.raw.vcf`\n";
 		print PARTSHELL "$GATK -T VariantFiltration -R $GENOME --filterExpression \"QD < 20.0 || ReadPosRankSum < -8.0 || FS > 10.0 || QUAL < \$MEANQUAL\" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --variant $PREDir/$k/$k.$j.concordance.raw.vcf --logging_level ERROR -o $PREDir/$k/$k.$j.concordance.flt.vcf\n";
 		print PARTSHELL "grep -v \"Filter\" $PREDir/$k/$k.$j.concordance.flt.vcf >$PREDir/$k/$k.$j.confidence.raw.vcf\n";
 		print PARTSHELL "#BQSR,ReduceReads\n";
 		print PARTSHELL "$GATK -T BaseRecalibrator -R $GENOME -I $PREDir/$k/$k.$j.realign.bam -knownSites $PREDir/$k/$k.$j.confidence.raw.vcf -o $PREDir/$k/$k.$j.recal_data.grp\n";
 		print PARTSHELL "$GATK -T PrintReads  -R $GENOME -I $PREDir/$k/$k.$j.realign.bam -BQSR $PREDir/$k/$k.$j.recal_data.grp -o $PREDir/$k/$k.$j.recal.bam\n";
 		print PARTSHELL "$GATK -T ReduceReads  -R $GENOME -I $PREDir/$k/$k.$j.recal.bam -o $PREDir/$k/$k.$j.reduced.bam\n";
 		print CLEAN "rm -rf $PREDir/$k/$k.$j.sort.bam*  $PREDir/$k/$k.$j.head.bam* $PREDir/$k/$k.$j.dedup.bam* $PREDir/$k/$k.$j.dedup.metrics* $PREDir/$k/$k.$j.realn* $PREDir/$k/$k.$j.realign* $PREDir/$k/$k.$j.gatk.raw.vcf* $PREDir/$k/$k.$j.samtools.raw.vcf* $PREDir/$k/$k.$j.concordance* $PREDir/$k/$k.$j.recal_data.grp* $PREDir/$k/$k.$j.recal.bam* $PREDir/$k/$k.$j.realign* $PREDir/$k/$k.$j.recal*\n";
 		
 		print QSUBPRE "sh $PARTSHELLDir/$k.$j.sh\n";
 		}
 	}
 
 
 		close BWA;
 		close Chrbam;
 		close PARTSHELL;
 	open QSUBSNPDB ,">$PREDir/qsubsnpdb.sh" or die $!;
 	open QSUBCALLSNP ,">$CALLSNPshellDir/qsubcallsnp.sh" or die $!;
 	open QSUBCALLINDEL,">$CALLINDELshellDir/qsubcallindel.sh" or die $!;
 	
 	foreach my $j(@chr){
 	open Bamlist ,">$CALLSNPLISTDir/$j.bam.list" 	or die $!;
 	open Chrsnpdb,">$PREDir/$j.snpdb.sh"	or die $!;
 	open Chrcallsnp,">$CALLSNPshellDir/$j.sh" or die $!;
 	open Chrcallindel,">$CALLINDELshellDir/$j.sh" or die $!;
 	
 	print QSUBSNPDB "sh $PREDir/$j.snpdb.sh\n";
 	print QSUBCALLSNP "sh $CALLSNPshellDir/$j.sh\n";
 	print QSUBCALLINDEL "sh $CALLINDELshellDir/$j.sh\n";
 	for my $k(sort keys %KEYNAMES){
 	print Bamlist "$PREDir/$k/$k.$j.reduced.bam\n";
}
	my $v=();
	for my $k(sort keys %KEYNAMES){
	$v .= "-V $PREDir/$k/$k.$j.confidence.raw.vcf ";
	}	
	
 	print Chrsnpdb "$GATK  -T CombineVariants -R $GENOME $v  -genotypeMergeOptions UNIQUIFY -o $PREDir/$j.snpdb.vcf\n";
 	print Chrcallsnp "$GATK -T UnifiedGenotyper -stand_call_conf 50.0 -stand_emit_conf 30  -rf BadCigar -glm BOTH -nt 2 -nct 4 -R $GENOME -I $CALLSNPLISTDir/$j.bam.list  --dbsnp $PREDir/$j.snpdb.vcf -o $CALLSNPDir/$j.raw.snps.indels.vcf\n";
 	print Chrcallsnp "$GATK -T VariantRecalibrator -R $GENOME -input $CALLSNPDir/$j.raw.snps.indels.vcf   -recalFile $CALLSNPDir/$j.output.recal -tranchesFile $CALLSNPDir/$j.output.tranches -mode SNP  -resource:concordantSet,VCF,known=true,training=true,truth=true,prior=10.0 $PREDir/$j.snpdb.vcf -an QD  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP\n";
  print Chrcallsnp "$GATK -T ApplyRecalibration  -R $GENOME -input $CALLSNPDir/$j.raw.snps.indels.vcf -tranchesFile $CALLSNPDir/$j.output.tranches -recalFile $CALLSNPDir/$j.output.recal  -mode SNP --ts_filter_level 99.0 -o $CALLSNPDir/$j.recalibrated.vcf\n";
  print Chrcallsnp "grep -E 'PASS|#' $CALLSNPDir/$j.recalibrated.vcf >$CALLSNPDir/$j.recalibrated.filtered.vcf\n";
  print Chrcallsnp "$GATK -T VariantFiltration   -R $GENOME  --filterExpression \"QD < 2.0 || MQ < 40.0 || ReadPosRankSum < -8.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5\" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR --variant $CALLSNPDir/$j.recalibrated.filtered.vcf -o $CALLSNPDir/$j.filtered.vcf\n"; 
  print Chrcallsnp "grep -E 'PASS|#' $CALLSNPDir/$j.filtered.vcf >$CALLSNPDir/$j.snp.final.vcf\n";
  print Chrcallsnp "perl $BinDir/snp.vcf_to_genotype.pl -input $CALLSNPDir/$j.snp.final.vcf -output $CALLSNPDir/$j.snp.genotype\n";
  
  print Chrcallindel "$GATK -T VariantRecalibrator -R $GENOME -input $CALLSNPDir/$j.raw.snps.indels.vcf --maxGaussians 3 -std 10.0 -percentBad 0.10 -recalFile $CALLINDELDir/$j.output.recal -tranchesFile $CALLINDELDir/$j.output.tranches -mode INDEL -resource:concordantSet,VCF,known=true,training=true,truth=true,prior=10.0 $PREDir/$j.snpdb.vcf -an QD  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP\n";
  print Chrcallindel "$GATK -T ApplyRecalibration  -R $GENOME -input $CALLSNPDir/$j.raw.snps.indels.vcf -recalFile $CALLINDELDir/$j.output.recal -tranchesFile $CALLINDELDir/$j.output.tranches -mode INDEL --ts_filter_level 95.0 -o $CALLINDELDir/$j.recalibrated.vcf\n";
	print Chrcallindel "grep -E 'PASS|#' $CALLINDELDir/$j.recalibrated.vcf >$CALLINDELDir/$j.recalibrated.filtered.vcf\n";
	print Chrcallindel "$GATK -T VariantFiltration   -R $GENOME --filterExpression \"QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0\" --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR --variant $CALLINDELDir/$j.recalibrated.filtered.vcf -o $CALLINDELDir/$j.filtered.vcf\n";
	print Chrcallindel "grep -E 'PASS|#' $CALLINDELDir/$j.filtered.vcf >$CALLINDELDir/$j.indel.final.vcf\n";
	print Chrcallindel "perl $BinDir/indel.vcf_to_genotype.pl -input $CALLINDELDir/$j.indel.final.vcf -output $CALLINDELDir/$j.indel.genotype\n";
	print CLEAN "rm -rf $PREDir/$j.snpdb.vcf* $CALLSNPDir/$j.recalibrated.filtered* $CALLSNPDir/$j.raw* $CALLSNPDir/$j.output.recal* $CALLSNPDir/$j.output.tranches* $CALLINDELDir/$j.output.recal* $CALLINDELDir/$j.output.tranches* $CALLSNPDir/$j.recalibrated.vcf* $CALLINDELDir/$j.recalibrated* $CALLSNPDir/$j.filtered.vcf* $CALLINDELDir/$j.filtered.vcf*\n"; 
	 }
	 close Bamlist;
	 close Chrsnpdb;
	 close PARTSHELL;
	 close CLEAN;
	

open  FINAL, ">$OUTDIR/final.sh" or die $!;
print FINAL "#!/bin/bash\n\n";
print FINAL "echo start at `date +'%Y-%m-%d %H:%M:%S %z'`\n\n";
print FINAL "# build index \n";
print FINAL "perl $Bin/qsub-sge.pl  --getmem   --lines 3 --convert no --resource vf=2G --queue $queue --subprjctid ngb_un  --jobprefix INDEX $INDEXDir/bulid_index.sh\n";
print FINAL "echo -e '\\t' build index done at `date +'%Y-%m-%d %H:%M:%S %z'`\n\n";
print FINAL "#BWA \n";
print FINAL "perl $Bin/qsub-sge.pl  --getmem  --maxjob 100 --convert no --reqsub --queue $queue  --subprjctid ngb_un  --resource vf=`du -b $GENOME | awk '{sum += \$1} END {printf(\"%.2fG\", sum / 1024 ^ 3*6)}'` --jobprefix BWA $BWADir/bwa.sh\n";
print FINAL "echo -e '\\t' bwa done at `date +'%Y-%m-%d %H:%M:%S %z'`\n\n";
print FINAL "#split bwa result \n";
print FINAL "perl $Bin/qsub-sge.pl  --getmem  --maxjob 100 --convert no --reqsub --queue $queue  --subprjctid ngb_un  --resource vf=`du -b $GENOME | awk '{sum += \$1} END {printf(\"%.2fG\", sum / 1024 ^ 3*10)}'` --jobprefix Chrbam $ChrbamDir/Chrbam.sh\n";
print FINAL "echo -e '\\t' split done at `date +'%Y-%m-%d %H:%M:%S %z'`\n\n";
print FINAL "#generate  preprocess result \n";
print FINAL "perl $Bin/qsub-sge.pl  --getmem  --maxjob 100 --convert no  --reqsub --queue $queue --subprjctid ngb_un  --resource vf=`du -b $GENOME | awk '{sum += \$1} END {printf(\"%.2fG\", sum / 1024 ^ 3*20)}'` --jobprefix PRE $PARTSHELLDir/qsubpre.sh\n";
print FINAL "echo -e '\\t' preprocess result done at `date +'%Y-%m-%d %H:%M:%S %z'`\n\n";
print FINAL "#generate  confidence snpdb file \n";
print FINAL "perl $Bin/qsub-sge.pl  --getmem  --maxjob 100  --convert no --reqsub --queue $queue --subprjctid ngb_un --resource vf=10G --jobprefix CON $PREDir/qsubsnpdb.sh\n";
print FINAL "echo -e '\\t' generate snpdb file done at `date +'%Y-%m-%d %H:%M:%S %z'`\n\n";
print FINAL "#call  snp and VQSR \n";
print FINAL "perl $Bin/qsub-sge.pl  --getmem  --maxjob 100 --convert no --reqsub  --queue $queue --subprjctid ngb_un --resource vf=`du -b $GENOME | awk '{sum += \$1} END {printf(\"%.2fG\", sum / 1024 ^ 3*20)}'` --jobprefix CALLSNP $CALLSNPshellDir/qsubcallsnp.sh\n";
print FINAL "echo -e '\\t'  generate snp result done at `date +'%Y-%m-%d %H:%M:%S %z'`\n\n";
print FINAL "#call  indel and VQSR \n";
print FINAL "perl $Bin/qsub-sge.pl  --getmem  --maxjob 100 --convert no --reqsub  --queue $queue --subprjctid ngb_un  --resource vf=`du -b $GENOME | awk '{sum += \$1} END {printf(\"%.2fG\", sum / 1024 ^ 3*20)}'` --jobprefix CALLINDEL $CALLINDELshellDir/qsubcallindel.sh\n";
print FINAL "echo -e '\\t'  generate indel result done at `date +'%Y-%m-%d %H:%M:%S %z'`\n\n";


close FINAL;
&showLog("done");

exit 0;
sub showLog {
	my ($info) = @_;
	my @times = localtime; # sec, min, hour, day, month, year
	print STDERR sprintf("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900, $times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}
sub MakeDir
{
	system("mkdir -m 755 -p $INDEXDir") if (!-d $INDEXDir);
	system("mkdir -m 755 -p $InFileDir") if (!-d $InFileDir);
	system("mkdir -m 755 -p $BWADir") if (!-d $BWADir);
	system("mkdir -m 755 -p $ChrbamDir") if (!-d $ChrbamDir);
	system("mkdir -m 755 -p $PREDir") if (!-d $PREDir);
	system("mkdir -m 755 -p $PARTSHELLDir") if (!-d $PARTSHELLDir);
	system("mkdir -m 755 -p $CALLSNPDir") if (!-d $CALLSNPDir);
	system("mkdir -m 755 -p $CALLSNPLISTDir") if (!-d $CALLSNPLISTDir);
	system("mkdir -m 755 -p $CALLSNPshellDir") if (!-d $CALLSNPshellDir);
	system("mkdir -m 755 -p $CALLINDELDir") if (!-d $CALLINDELDir);
	system("mkdir -m 755 -p $CALLINDELshellDir") if (!-d $CALLINDELshellDir);
}
	 
	                              
	
	
	
	
