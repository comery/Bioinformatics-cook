#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
use lib $Bin;
use Support_Program;
use Cwd 'abs_path';

sub usage {
	print <<USAGE;
usage:
	perl $0 [options]
author
	nixiaoming   nixiaoming\@genomics.cn
description:

options:
	-help        : print help info
	-lib    (str): 
	-outdir (str): (default = . )
	-queue1  (str): (default = bc_gag.q,all.q )
	-queue2  (str): (default = bc_gag.q,all.q  )  # with 32G Compute Node
	-fqdir  (str): (default = ./fq )
	-path1  (str): (default = $Bin )
	-path2  (str): (default = /ifs1/EPI/tangqq/bin/common/ )
e.g.:
	perl $0 -lib denovo.lib 
	perl $0 -lib denovo.lib -outdir denovo -fqdir fq  -queue1 bc_gag.q,all.q -queue2 bc_gag.q,all.q  -path1 $Bin -path2 /ifs1/EPI/tangqq/bin/common 
USAGE
}
##########################################################################################
############ 初始化，变量
##########################################################################################
my $path1=$Bin;
my $rootpath='/ifs4/BC_PUB/biosoft/pipe/bc_ba'; #add by tangqq
my $rna_lib='/ifs4/BC_PUB/biosoft/pipe/bc_ba/software/RNA_lib/lib';
my $path2 ="$rootpath/software";
my $queue1 ='bc.q,all.q';
my $queue2 ='bc.q,all.q';
my ($help,$lib,$outdir,$fqdir);
GetOptions(
	"help"=>\$help,
	"lib=s"=>\$lib,
	"outdir=s"=>\$outdir,
	"queue1=s"=>\$queue1,
	"queue2=s"=>\$queue2,
	"fqdir=s"=>\$fqdir,
	"path1=s"=>\$path1,
	"path2=s"=>\$path2,
);


if ((defined $help) || (! defined $lib) ){
	&usage();
	exit 0;
}
$outdir ||='.';
system "mkdir -p $outdir";
$outdir=abs_path($outdir);
$fqdir ||='./fq';
system "mkdir -p $fqdir";
$fqdir=abs_path($fqdir);
$lib = abs_path($lib);


#########################################################################
my $shdir=$outdir.'/SH';
system "mkdir -p $shdir" unless(-d $shdir);
my $samsh=$shdir.'/SampleSH';
system "mkdir -p $samsh" unless(-d $samsh);
my $oldshdir=$shdir.'/oldsh';
system "mkdir -p $oldshdir " unless(-d $oldshdir);
my $cwdmv='find  '.$shdir.' -maxdepth 1 -name "*.sh.*"  -exec mv {} '.$oldshdir.' ";"';
system "$cwdmv";
my $version_dir=$outdir.'/version';
######################added by dushuo 2010-11-15
my $versionsh=$outdir.'/version.sh';
open VERSION,'>',$versionsh;
print VERSION 'export LD_LIBRARY_PATH="'.$rna_lib.':$LD_LIBRARY_PATH"';
print VERSION "\n";
#############marked by dushuo
system "mkdir  -p $version_dir" unless (-d $version_dir);
my %config;
my $pro_log="$version_dir\/program_path.txt";
Support_Program(\%config,$path1,$path2,1,$pro_log); #1表示 debug 
my $cwd='perl  '.$config{"svn_version"}.'  '.$Bin.'  '.$0.'     '.$config{"svn_version"}.'  > '.$version_dir.'/denovo.txt '."\n";
system "$cwd";
my $time='`date +%y-%m-%d.%H:%M:%S`';

my $for_clu_file='';
my $for_qua_file='';###########added by dushuo
my $libdir=$outdir.'/lib';
system "mkdir -p $libdir" unless(-d $libdir);
my $unigene_file='';
my $upload_sh=$shdir.'/upload.sh';
open UPLOADSH,'>',$upload_sh or die "can't open the upload sh $upload_sh";
my $clear_sh=$shdir.'/clear.sh';
open CLEARSH,'>',$clear_sh or die "can't open the  sh for rm the the mid file $clear_sh ";


my $cat_qua_contig='';
my $cat_qua_scaf='';
my $cat_qua_f1='';
##my $cat_qua_f2='';
my $cat_qua_fq='';

################ 读入lib
my %option=();
my @samplename=();
open IN,$lib or die "can't open the file $lib";
$/='>';<IN>;$/="\n";
while (my $sample=<IN>) {
	chomp $sample;
	$sample=~s/\s*$//;
	push (@samplename,$sample);
	$/="\n>";
	my $block=<IN>;
	chomp $block;
	$/="\n";
	my @Lines=split /\n/,$block;
	foreach my $line (@Lines) {
		next if($line=~/^\s*$/);
		next if ($line=~/^\s*\#/);
		$line=~s/\s*$//;
		if ($line=~/^(\S+)\s*=\s*(\S.*)$/) {
			$option{$sample}{$1}=$2;
		}
		if ($line=~/^q[12]?\s*=\s*\S+$/) {
			$option{$sample}{fq}.=$line."\n";
		}
		if ($line=~/^length\s*=\s*(\d+)\,\s*(\d+)\s*$/){
			$option{$sample}{length} = ($1<$2)?$2:$1;
		}
	}
}
close IN;
#################################################
my $cluster=shift @samplename;
my $sam_num=@samplename;###########added by dushuo for check.pl
my $upload_dir=$outdir.'/'.$cluster;
system "mkdir -p $upload_dir" unless(-d $upload_dir);

#############################
# warnings for project_name && sub_project_name
# 检测是否输入 项目名称等
############################
#unless(exists $option{$cluster}{'denovo_info'}){
#        die 'Please specify the file path, which used to collect the project informations! It is very importent!',"\n";
#}
unless(exists $option{$cluster}{'pro_name'}){
        die 'Please input the project\'s name, which you can get it from the BMS web-page! For example: pro_name=xxxx',"\n";
}
unless(exists $option{$cluster}{'sub_prjct_id'}){
	$option{$cluster}{'sub_prjct_id'} = "batest";
}
foreach my $sssam(@samplename){
        unless(exists $option{$sssam}{'name'}){
                die 'Please input the sub_project\'s name, which you can get it from the BMS web-page! For example: name=xxx. For your sample: ',$sssam,"\n";
        }
}

my $up_as_fi_dir=$upload_dir.'/assembly';
my $up_contig_dir=$up_as_fi_dir.'/1.Contig';
#my $up_scaf_dir=$up_as_fi_dir.'/2.Scaffold';
my $up_clu_dir=$up_as_fi_dir.'/2.Unigene';
##########################################################################################
#############  组装 聚类
##########################################################################################
######################## 打开句柄

foreach my $sam (@samplename) {
######################################################### 生成 lib
	my $samplelib=$libdir.'/'.$sam.'.lib';
	open LIB,'>',$samplelib or die "can't open the $sam lib $samplelib";
	if (! exists $option{$sam}{name} ) {
		$option{$sam}{name}=$sam;
	}
	print LIB 'name='.$option{$sam}{name}."\n";
	if (exists $option{$sam}{reverse_seq}) {
		print LIB 'reverse_seq='.$option{$sam}{reverse_seq}."\n";
	}else{
		print LIB 'reverse_seq=0'."\n";
	}
	if (exists $option{$sam}{rank}) {
		print LIB 'rank='.$option{$sam}{rank}."\n";
	}
##################################added by wanshun at 20110811###############################################
	if (exists $option{$sam}{'jaccard_clip'}) {
		print LIB 'jaccard_clip='.$option{$sam}{'jaccard_clip'}."\n";
	}
	if (exists $option{$sam}{'bfly_opts'}) {
		print LIB 'bfly_opts='.$option{$sam}{'bfly_opts'}."\n";
	}
	if (exists $option{$sam}{'seqType'}) {
		print LIB 'seqType='.$option{$sam}{'seqType'}."\n";
	}
	if (exists $option{$sam}{'min_contig_length'}) {
		print LIB 'min_contig_length='.$option{$sam}{'min_contig_length'}."\n";
	}
	if (exists $option{$sam}{'CPU'}) {
		print LIB 'CPU='.$option{$sam}{'CPU'}."\n";
	}
	if (exists $option{$sam}{'min_kmer_cov'}) {
                print LIB 'min_kmer_cov='.$option{$sam}{'min_kmer_cov'}."\n";
        }
	if (exists $option{$sam}{'SS_lib_type'}) {
		print LIB 'SS_lib_type='.$option{$sam}{'SS_lib_type'}."\n";
	}
	if (exists $option{$sam}{'path_reinforcement_distance'}) {
		print LIB 'path_reinforcement_distance='.$option{$sam}{'path_reinforcement_distance'}."\n";
	}
	if (exists $option{$sam}{'group_pairs_distance'}) {
		print LIB 'group_pairs_distance='.$option{$sam}{'group_pairs_distance'}."\n";
	}
	if (exists $option{$sam}{'min_glue'}) {
		print LIB 'min_glue='.$option{$sam}{'min_glue'}."\n";
	}
##############################生成\$sam.sh,用于并行处理####并行相关changed by dushuo 不逐一注释
	my $samshdir="$samsh/$sam";
	system "mkdir -p  $samshdir " unless(-d $samshdir);
	my $sam_sh=$samsh."\/$sam.sample_assembly.sh";
	open SAMSH,'>',$sam_sh;
##############################处理fq 
	if (defined $option{$sam}{fq}) {
		my $tmpfqdir="$outdir/fq/$sam";
		system "mkdir -p $tmpfqdir/nodup" unless (-d "$tmpfqdir/nodup");
		my @fq=split /\n/,$option{$sam}{fq};
		$option{$sam}{fq}='';
		my $num_fq=1;
		my $i=0;
		my %fqindir=();###存放原始fq的路径；
		my $libfq='';
		my $filterfq='';
		while ($i<@fq) {
			chomp $fq[$i];
			my $tmpfq="$tmpfqdir/$sam\_l$num_fq";
			if ($fq[$i]=~/^q[12]\s*=\s*\S+_[12]\.(?:fq\.gz|fq)/) {#####找出成对的fq 或 fq.gz 
				my ($fq1,$fq2);
				if ($fq[$i]=~/^q1\s*=\s*(\S+_1\.)(fq\.gz|fq)\s*$/) {
					$fq1=$1.$2;
					chomp $fq[$i+1];
					if($fq[$i+1]=~/^q2\s*=\s*(\S+_2\.)(fq\.gz|fq)\s*$/){
						$fq2=$1.$2;
					}
					$i++;
				}elsif($fq[$i]=~/^q2\s*=\s*(\S+_2\.)(fq\.gz|fq)\s*$/){
					$fq1=$1.$2;
					chomp $fq[$i+1];
					if($fq[$i+1]=~/^q1\s*=\s*(\S+_2\.)(fq\.gz|fq)\s*$/){
						$fq2=$1.$2;
					}
					my $aaa=$fq2;
					$fq2=$fq1;
					$fq1=$aaa;
					$i++;
				}else{
					print STDERR "warnings : error in lib of $sam fq \n $fq[$i]\n $fq[$i+1]\n";
				}
				$i++;
				if (exists $option{$sam}{"filter_fq_option"}) {
					my $sam_filter="$samshdir/$sam.filter.sh";
					open SAMFILTER,'>',$sam_filter;
					print SAMFILTER "export LD_LIBRARY\_PATH=$rna_lib:\$LD_LIBRARY_PATH\n";
					print SAMFILTER $config{"filter_fq"}.' --fq1 '.$fq1.'  --fq2 '.$fq2.'   --out '.$tmpfq.'    '.$option{$sam}{"filter_fq_option"}.'  && echo finish filter '.$tmpfq.' at time '.$time."\n";
					my $suffix="fq".fq_suffix($sam,$option{$sam}{"filter_fq_option"});
					print SAMFILTER "mv $tmpfq\_1.$suffix  $tmpfq\_2.$suffix  $fqdir \n";
					print SAMFILTER "ln -s  $fqdir/$sam\_l$num_fq\_1.$suffix  $fqdir/$sam\_l$num_fq\_2.$suffix    $tmpfqdir/  \n";
					print SAMFILTER "ln -s  $tmpfq.stat  $tmpfqdir/nodup/$sam\_l$num_fq.stat    \n";
					print SAMFILTER "cd $fqdir\n";
					print SAMFILTER $config{"md5sum"}." $sam\_l$num_fq\_1.$suffix  > $sam\_l$num_fq\_1.$suffix.md5  \n";
					print SAMFILTER $config{"md5sum"}." $sam\_l$num_fq\_2.$suffix  > $sam\_l$num_fq\_2.$suffix.md5  \n\n";
					print VERSION 'perl  '.$config{svn_version}.'     '.$config{"filter_fq"}.'   '.$config{qsub_sge}.'  > '.$version_dir.'/filter_fq.txt '."\n";
					$filterfq.="q1=$tmpfq\_1.$suffix\n";
					$filterfq.="q2=$tmpfq\_2.$suffix\n";
					#####################################################added by dushuo
					if (exists $option{$sam}{dup_fq_option}){
						my $rate.=$option{$sam}{dup_fq_option};
						print SAMFILTER "perl  ".$config{"duplication.gz.pl"}." -fq1 $tmpfq\_1.$suffix -fq2 $tmpfq\_2.$suffix -out $tmpfqdir/nodup/$sam\_l$num_fq  $rate \n";
						print VERSION 'perl  '.$config{svn_version}.'  '.'    '.$config{"duplication.gz.pl"}.'  '.$config{"duplication_single.gz.pl"}.'  '.$config{qsub_sge}.'  > '.$version_dir.'/dup_fq.txt '."\n";
					}else{
						print SAMFILTER "perl  ".$config{"duplication.gz.pl"}." -fq1 $tmpfq\_1.$suffix -fq2 $tmpfq\_2.$suffix -out $tmpfqdir/nodup/$sam\_l$num_fq \n";
						print VERSION 'perl  '.$config{svn_version}.'  '.'    '.$config{"duplication.gz.pl"}.'  '.$config{"duplication_single.gz.pl"}.'  '.$config{qsub_sge}.'  > '.$version_dir.'/dup_fq.txt '."\n";
					}
#					print SAMSH "sh $sam_filter\n\n";
#					print SAMSH<<SAMCHECK;
#if [ \$? -ne 0 ]
#then
#	echo error in $sam.sample_assembly -- filter
#	echo error in $sam.sample_assembly -- filter 1>&2
#	exit 1
#fi
#
#SAMCHECK
				}else{
					my $sam_filter="$samshdir/$sam.filter.sh";
					open SAMFILTER,'>',$sam_filter;
					print SAMFILTER "perl  ".$config{"duplication.gz.pl"}." -fq1 $fq1 -fq2 $fq2 -out $tmpfqdir/nodup/$sam\_l$num_fq \n";
					print VERSION 'perl  '.$config{svn_version}.'  '.'    '.$config{"duplication.gz.pl"}.'  '.$config{"duplication_single.gz.pl"}.'  '.$config{qsub_sge}.'  > '.$version_dir.'/dup_fq.txt '."\n";
					$filterfq.="q1=$fq1\n";
					$filterfq.="q2=$fq2\n";
					close SAMFILTER;
					print SAMSH "sh $sam_filter\n\n";
					print SAMSH<<SAMCHECK;
if [ \$? -ne 0 ]
then
	echo error in $sam.sample_assembly -- filter
	echo error in $sam.sample_assembly -- filter 1>&2
	exit 1
fi

SAMCHECK
				}
				$libfq.="q1=$tmpfqdir/nodup/$sam\_l$num_fq\_1.fq\n";
				$libfq.="q2=$tmpfqdir/nodup/$sam\_l$num_fq\_2.fq\n";
				print CLEARSH "rm    $tmpfqdir/nodup/$sam\_l$num_fq\_1.fq   $tmpfqdir/nodup/$sam\_l$num_fq\_2.fq $tmpfqdir/nodup/$sam\_l$num_fq.dup.list  && echo finish rm $tmpfq\n";
			}elsif($fq[$i]=~/^q\s*=\s*(\S+\.)(fq\.gz|fq)\s*$/){#####单个的fq
				my $fq=$1.$2;
				$i++;
				my $fqname=basename($fq);
				####过滤单个的fq
				if (exists $option{$sam}{filter_fq_option}) {
					my $sam_filter="$samshdir/$sam.filter.sh";
					open SAMFILTER,'>',$sam_filter;
					print SAMFILTER $config{"filter_fq"}.' --fq1 '.$fq.' --out '.$tmpfq.'    '.$option{$sam}{"filter_fq_option"}.'  && echo finish filter '.$tmpfq.' at time '.$time."\n";
					my $suffix="fq".fq_suffix($sam,$option{$sam}{"filter_fq_option"});
					print SAMFILTER "mv $tmpfq.$suffix   $fqdir \n";
					print SAMFILTER " ln -s  $fqdir/$sam\_l$num_fq.$suffix   $tmpfqdir/$sam\_l$num_fq.$suffix  \n";
					print SAMFILTER " ln -s  $tmpfq.stat    $tmpfqdir/nodup/$sam\_l$num_fq.stat  \n";
					print SAMFILTER "cd $fqdir\n";
					print SAMFILTER $config{"md5sum"}." $sam\_l$num_fq.$suffix  > $sam\_l$num_fq.$suffix.md5  \n";
					print SAMFILTER "echo \n\n";
					print SAMFILTER "perl  ".$config{"duplication_single.gz.pl"}." -fq $tmpfq.$suffix   -out $tmpfqdir/nodup/$sam\_l$num_fq \n";
					
					  24     -queue2  (str): (default = bc_gag.q,all.q  )  # wprint VERSION 'perl  '.$config{svn_version}.'     '.$config{"filter_fq"}.'   '.$config{qsub_sge}.'  > '.$version_dir.'/filter_fq.txt '."\n";
					print VERSION 'perl  '.$config{svn_version}.'  '.'    '.$config{"duplication.gz.pl"}.'  '.$config{"duplication_single.gz.pl"}.'  '.$config{qsub_sge}.'  > '.$version_dir.'/dup_fq.txt '."\n";
					$filterfq.="q=$tmpfq.$suffix\n";
					close SAMFILTER;
					print SAMSH "sh $sam_filter\n";
					print SAMSH<<SAMCHECK;
if [ \$? -ne 0 ]
then
	echo error in $sam.sample_assembly
	echo error in $sam.sample_assembly 1>&2
	exit 1
fi

SAMCHECK
				}else{
					my $sam_filter="$samshdir/$sam.filter.sh";
					open SAMFILTER,'>',$sam_filter;
					print SAMSH "perl  ".$config{"duplication_single.gz.pl"}." -fq $fq   -out $tmpfqdir/nodup/$sam\_l$num_fq \n";
					print VERSION 'perl  '.$config{svn_version}.'  '.'    '.$config{"duplication.gz.pl"}.'  '.$config{"duplication_single.gz.pl"}.'  '.$config{qsub_sge}.'  > '.$version_dir.'/dup_fq.txt '."\n";
					$filterfq.="q=$fq\n";
					close SAMFILTER;
					print SAMSH "sh $sam_filter\n\n";
					print SAMSH<<SAMCHECK;
if [ \$? -ne 0 ]
then
	echo error in $sam.sample_assembly -- filter
	echo error in $sam.sample_assembly -- filter 1>&2
	exit 1
fi

SAMCHECK
				}

				####fq 写入lib
				$libfq.="q=$tmpfqdir/nodup/$sam\_l$num_fq.fq\n";
				print CLEARSH 'rm  '."  $tmpfqdir/nodup/$sam\_l$num_fq.dup.list  && echo rm $tmpfqdir/nodup/$sam\_l$num_fq.dup.list\n";
			}else{
				die "error in lib of $sam fq ";
			}
			$num_fq++;
		}
		print LIB $libfq;
		$option{$sam}{fq}=$filterfq;
	}else{
		die "please input  $sam  fq file \n";
	}	
	close LIB;
	$cat_qua_fq.='   '.$samplelib;
######################################################### 组装 
	my $grapedir=$outdir.'/assembly_fill';
	system "mkdir -p   $grapedir" unless(-d $grapedir);
	my $sample_dir=$grapedir.'/'.$sam;
	system "mkdir -p  $sample_dir" unless(-d $sample_dir);
	my $sample_as_dir=$sample_dir.'/assembly';
	system "mkdir -p  $sample_as_dir" unless(-d $sample_as_dir);
	my ($contigname,$scafname);
	if (!exists $option{$sam}{"map_len"}){
		$option{$sam}{"map_len"} = 50;
	}
	if (!exists $option{$sam}{"contig_len"}) {
		$option{$sam}{"contig_len"} = 100;
	}
	if (!exists $option{$sam}{"scaffold_len"}) {
		$option{$sam}{"scaffold_len"} = 100;
	}
	if (!exists $option{$sam}{"grape_option"}){
		$option{$sam}{"grape_option"} = "-K 25  -M 2  -D 1";
	}
	if (!exists $option{$sam}{"fill_option"}){
		$option{$sam}{"fill_option"} = " -p 27  -t 6";
	}

	my $sam_as="$samshdir/$sam.as.sh";
	if((exists  $option{$sam}{contig})){
		$option{$sam}{contig}=`perl $config{absolute_dir} $option{$sam}{contig}`;
		unless(-s $option{$sam}{contig}) {
			die 'warning  : '.$option{$sam}{contig}.' is not exists or it size is zero '."\n";
		}
		open SAMAS,'>',$sam_as;
		print SAMAS 'export LD_LIBRARY_PATH="'.$rna_lib.':$LD_LIBRARY_PATH"';
		print SAMAS "\n";
		$contigname=basename($option{$sam}{contig});
	}
	elsif ((exists $option{$sam}{"bfly_opts"}) || (exists $option{$sam}{"jaccard_clip"})) {
		my $sam_as="$samshdir/$sam.as.sh";
		open SAMAS,'>',$sam_as;
		print SAMAS 'export LD_LIBRARY_PATH="'.$rna_lib.':$LD_LIBRARY_PATH"';
                print SAMAS "\n";
		my $K=31;
#		print SAMAS "sh ".$config{"Trinity"}." -s $samplelib -o $sample_as_dir/$sam -n $sam".' &&  echo  finish '.$sam.' grape  at time '.$time."\n"; #add by tangqq 2011.6.9 21.07
		print SAMAS<<CHECKTR;
if [ ! -f $sample_as_dir/$sam/Trinity.fasta ];then
	echo No result of Trinity!
	exit 1
fi
CHECKTR
		print VERSION 'perl  '.$config{svn_version}.'  '.$config{get_chosen_fa}.'  '.$config{fa_quality}.'  '.$config{barplot}.'  '.$config{qsub_sge}.'  > '.$version_dir.'/assemble.txt '."\n";#by chuanjun
		$contigname=$sam.'.contig';
		$scafname=$sam.'.scafSeq';
		print CLEARSH 'rm -rf '.$sample_as_dir.'/'.$sam.'/'.'both.fa* '.$sample_as_dir.'/'.$sam.'/'.'left.fa '.$sample_as_dir.'/'.$sam.'/'.'right.fa '.$sample_as_dir.'/'.$sam.'/'.'chrysalis '.$sample_as_dir.'/'.$sam.'/'.'jaccard_clip_workdir* '.$sample_as_dir.'meryl* '.$sample_as_dir.'/'.$sam.'/'.'mer_counts_0'.'  && echo finish rm the mid assembly file '.$sam.' at time '.$time."\n";
	}else{
			print STDERR 'warnings : no assembling  of '.$sam."\n".' you can add grape_option=xxxx in lib file to  assembling or give the  assembled fa file in lib file by add contig= xxxx \n scaffold= xxx'."\n";
	}
	if((exists  $option{$sam}{contig}) || (exists $option{$sam}{grape_option})){
		print SAMAS 'perl '.$config{get_chosen_fa}.' -fa '.$sample_as_dir.'/'.$contigname.'   -output '.$sample_as_dir.'/'.$sam.'-Contig.fa   -type Contig -name '.$sam.' -size  -len '.$option{$sam}{"contig_len"}.' && echo finish get the length contig and changed the contig  names'."\n";  #remove by tangqq at 2011.6.9 21:16
		print SAMAS "ln -s $sample_as_dir/$sam/$sam-Contig.fa $sample_as_dir\n";
		print SAMAS "ln -s $sample_as_dir/$sam/$sam.GapFilling.fa $sample_as_dir\n";
		print SAMAS 'perl '.$config{fa_quality}.' -len -Head -gap -N  -gc '.$sample_as_dir.'/'.$sam.'-Contig.fa    &&  echo finish statistics  contig  scaffold for draw pdf at time '.$time."\n";
		print SAMAS ' perl '.$config{barplot}.'   '.$sample_as_dir.'/'.$sam.'-Contig.fa.quality.xls   '.$sam.'-Contig   && echo finish draw  pdf of contig at time '.$time."\n";
		print SAMAS "\n\n";
		close SAMAS;
		print SAMSH "sh $sam_as\n\n";
		print SAMSH<<SAMCHECK;
if [ \$? -ne 0 ]
then
	echo error in $sam.sample_assembly -- sam_as
	echo error in $sam.sample_assembly -- sam_as 1>&2
	exit 1
fi

SAMCHECK

		system "mkdir -p  $up_as_fi_dir" unless(-d $up_as_fi_dir);
		system "mkdir -p  $up_contig_dir " unless(-d $up_contig_dir);
		my @as_contig;
		push @as_contig,"$sample_as_dir/$sam-Contig.length.svg","$sample_as_dir/$sam-Contig.fa","$sample_as_dir/$sam-Contig.length.png","$sample_as_dir/$sam-Contig.length.txt";
		foreach my $assem_result(@as_contig){
			my $assem_filename = basename($assem_result);
			my $assem_target = $up_contig_dir.'/'.$assem_filename;
			my $assem_shell = &_Upload($assem_result,$assem_target);
			print UPLOADSH $assem_shell;
		}
		undef @as_contig;
		$cat_qua_contig.=$sample_as_dir.'/'.$sam.'-Contig.fa.quality.xls   ';###############changed by dushuo
	}
  ####################################################### 补洞
	my $sample_fi_dir=$sample_dir.'/fillgap';
	system "mkdir -p  $sample_fi_dir" unless(-d $sample_fi_dir);
	my $f1name;
	if(exists $option{$sam}{f1_fa}) {
		my $sam_f1="$samshdir/$sam.f1.sh";
		open SAMF1,'>',$sam_f1;
		print SAMF1 'export LD_LIBRARY_PATH="'.$rna_lib.':$LD_LIBRARY_PATH"';
		print SAMF1 "\n";
		$option{$sam}{f1_fa}=`perl $config{absolute_dir} $option{$sam}{f1_fa}`;
		unless(-s $option{$sam}{f1_fa}) {
			die 'warning  : '.$option{$sam}{f1_fa}.' is not exists or it size is zero '."\n";
		}
		print SAMF1 'ln -s '.$option{$sam}{f1_fa}.'  '.$sample_fi_dir.'  && echo finish ln '.$sample_fi_dir.'  at time '.$time."\n";
		$f1name=basename($option{$sam}{f1_fa});
		$f1name=$sample_fi_dir.'/'.$f1name;
		$cat_qua_f1.=$f1name.'.quality.xls  ';#########changed by dushuo
		print SAMF1 'perl '.$config{fa_quality}.' -len -Head -gap -N  -gc '.$f1name.'     && echo finish  statistics quality  for draw pdf of '.$sam.'.GapFilling.fa  at time '.$time."\n";
		my $f1pdf=basename($f1name);
		$f1pdf=~s/\.fa$//;
		print SAMF1 ' perl '.$config{barplot}.'   '.$f1name.'.quality.xls   '.$f1pdf.'  -gap     && echo finish  draw pdf of '.$f1name.' at time '.$time."\n";
		print SAMF1 "\n";
	}elsif (exists $option{$sam}{fill_option}) {
		my $sam_f1="$samshdir/$sam.f1.sh";
		open SAMF1,'>',$sam_f1;
		print SAMF1 'export LD_LIBRARY_PATH="'.$rna_lib.':$LD_LIBRARY_PATH"';
        print SAMF1 "\n";
		print VERSION 'perl  '.$config{svn_version}.'      '.$config{fa_quality}.'  '.$config{barplot}.'  '.$config{qsub_sge}.'  > '.$version_dir.'/fillgap.txt '."\n"; # dushuo
		$f1name=$sample_fi_dir.'/'.$sam.'.GapFilling.fa';
		print SAMF1 "ln -s $sample_as_dir/$sam/$sam.GapFilling.fa  $sample_fi_dir\n";
		$cat_qua_f1.=$f1name.'.quality.xls  ';#########changed by dushuo
		my $f1pdf=basename($f1name);
		$f1pdf=~s/\.fa$//;
		print SAMF1 "\n";
		close SAMF1;
		print SAMSH "sh $sam_f1\n\n";
		print SAMSH<<SAMCHECK;
if [ \$? -ne 0 ]
then
	echo error in $sam.sample_assembly -- sam_f1
	echo error in $sam.sample_assembly -- sam_f1 1>&2
	exit 1
fi

SAMCHECK
	}else{
		print STDERR 'warnings : no fill gap   of '.$sam."\n".'you can add fill_option=xxxx in lib file to  fill gap or give the  filled  fa file in lib file by add f1_fa= xxxx'."\n";
	}

####################################################### sample 聚类
	my $sample_clu_dir=$sample_dir.'/cluster';
	my ($samgene,$samgenepdf);
	if (exists $option{$sam}{unigene}) {
		my $sam_clust="$samshdir/$sam.cluster.sh";
		open SAMCLU,'>',$sam_clust;
		print SAMCLU 'export LD_LIBRARY_PATH="'.$rna_lib.':$LD_LIBRARY_PATH"';
		print SAMCLU "\n";

		system "mkdir -p  $up_clu_dir " unless(-d $up_clu_dir);
		system "mkdir -p  $sample_clu_dir " unless(-d $sample_clu_dir);
		$option{$sam}{unigene}=`perl $config{absolute_dir} $option{$sam}{unigene}`;
		unless(-s $option{$sam}{unigene}) {
			die 'warning  : '.$option{$sam}{unigene}.' is not exists or it size is zero '."\n";
		}
		print SAMCLU 'ln -s '.$option{$sam}{unigene}.'  '.$sample_clu_dir.'  && echo finish ln '.$sample_as_dir.'  at time '.$time."\n";
		print SAMCLU 'echo '."\n";
		print SAMCLU 'echo '."\n";
		print SAMCLU 'echo '."\n";
		print SAMCLU 'echo '."\n";
		print SAMCLU 'echo '."\n";
		print SAMCLU 'echo '."\n";
		print SAMCLU 'echo '."\n";
		print SAMCLU 'echo '."\n";
		print SAMCLU 'echo '."\n";
		$samgene=basename($option{$sam}{unigene});
		$samgenepdf=$samgene;
		$samgenepdf=~s/\.fa$//;
		my @sam_clus;
		push @sam_clus,"$sample_clu_dir/$samgene","$sample_clu_dir/$samgenepdf.length.svg","$sample_clu_dir/$samgenepdf.length.png","$sample_clu_dir/$samgenepdf.length.txt";
		foreach my $sam_clu_result(@sam_clus){
			my $sam_filename = basename($sam_clu_result);
			my $sam_target = $up_clu_dir.'/'.$sam_filename;
			my $sam_shell = &_Upload($sam_clu_result,$sam_target);
			print UPLOADSH $sam_shell;
		}
		undef @sam_clus;
		$for_clu_file.=$sample_clu_dir.'/'.$samgene.'  ';
		$for_qua_file.=$sample_clu_dir.'/'.$samgene.'.quality.xls ';#########added by dushuo
		print SAMCLU 'perl '.$config{fa_quality}.' -len -Head -gap -N  -gc '.$sample_clu_dir.'/'.$samgene.'   && echo finish statistics  quality for draw pdf of '.$sam.'-Unigene.fa  at time '.$time."\n";
		print SAMCLU ' perl '.$config{barplot}.'   '.$sample_clu_dir.'/'.$samgene.'.quality.xls   '.$samgenepdf.'  -gap   && echo finish draw pdf of '.$sam.'-Unigene.fa at time '.$time."\n\n";
	}elsif (exists $option{$sam}{tgicl_option}) {
		my $sam_clust="$samshdir/$sam.cluster.sh";
		open SAMCLU,'>',$sam_clust;
		print SAMCLU 'export LD_LIBRARY_PATH="'.$rna_lib.':$LD_LIBRARY_PATH"';
		print SAMCLU "\n";
		system "mkdir -p  $sample_clu_dir " unless (-d $sample_clu_dir );
		$samgene=$sam.'-Unigene.fa';
		$samgenepdf=$samgene;
		$samgenepdf=~s/\.fa$//;
		#if ((exists $option{$sam}{cluster}) && ($option{$sam}{cluster} eq 'f2')) {
			#print SAMSH 'ln -s  '.$f2name.'     '.$sample_clu_dir.'/for_'.$sam.'.cluster.fa'."\n";
		#}else{
		print SAMCLU 'ln -s  '.$f1name.'     '.$sample_clu_dir.'/for_'.$sam.'.cluster.fa'."\n";
		#}
		$option{$sam}{"unigene_len"} ||=200;
		print SAMCLU <<TGICLSAM;
cd $sample_clu_dir && echo cd $sample_clu_dir
perl $config{"tgicl"} for_$sam.cluster.fa $option{$sam}{"tgicl_option"}  && echo finish tgicl at $time 
cat asm_*/align > align   && echo get align at $time  
perl $config{"phrap.id.list.pl"} align align && echo get cluster id  at $time  
perl $config{"get_single.pl"}  align.cluster for_$sam.cluster.fa single  && echo get single  at $time 
cat asm_*/contigs > asm_cluster.fa  && echo get cluster sequences at $time
cat single.prefect.fa >> single.fa  && echo get singleton sequences at $time
cat asm_*/contigs  single.fa > tgicl_cluster_and_single.fa && echo get cluster out file  at $time
cat align.cluster single.list > tgicl_cluster_and_single.fa.list && echo get tgicled unigenes list file at $time
/opt/blc/genome/bin/formatdb -p F -i tgicl_cluster_and_single.fa && echo formatdb tgicl_cluster_and_single.fa at $time
/opt/blc/genome/bin/blastall -p blastn -m 8 -e 1e-10 -F F -a 10 -d $sample_clu_dir/tgicl_cluster_and_single.fa -i $sample_clu_dir/tgicl_cluster_and_single.fa -o $sample_clu_dir/all_vs_all.blast.m8 && echo blastn all to all at $time
perl $config{"cluster_for_coverage.pl"} $sample_clu_dir/all_vs_all.blast.m8 $sample_clu_dir && echo cluster based on coverage at $time
perl $config{"clusterOrSingle.pl"} $sample_clu_dir $sam $option{$sam}{"unigene_len"} && echo rename cluster and chose multiple cluster at $time
perl $config{"get_chosen_fa"} -fa $sam.all.fa -output $sam-Unigene.fa -len $option{$sam}{"unigene_len"} && echo get $option{$sam}{"unigene_len"} unigene fa at $time
perl $config{"get_chosen_fa"} -fa cluster.all.fa -output $sam-Cluster.fa -len $option{$sam}{"unigene_len"} && echo get $option{$sam}{"unigene_len"} cluster fa at $time
perl $config{"get_chosen_fa"} -fa singleton.all.fa -output $sam-Single.fa -len $option{$sam}{"unigene_len"} && echo get $option{$sam}{"unigene_len"} singleton fa at $time
mkdir -p id_tmp
ln -s $sample_clu_dir/cluster_and_single.fa.id.list id_tmp/cluster_and_single.fa.id.list
perl $config{"fa_quality"} -len -Head -gap -N  -gc   $samgene && echo fa quality   at $time 
perl $config{"barplot"}  $sam-Unigene.fa.quality.xls   $samgenepdf  -gap  && echo plot   at $time
TGICLSAM
		close SAMCLU;
		print SAMSH "sh $sam_clust\n\n";
		print SAMSH<<SAMCHECK;
if [ \$? -ne 0 ]
then
	echo error in $sam.sample_assembly -- sam_clu
	echo error in $sam.sample_assembly -- sam_clu 1>&2
	exit 1
fi

SAMCHECK
		$for_clu_file.=$sample_clu_dir.'/'.$samgene.'  ';
		$for_qua_file.=$sample_clu_dir.'/'.$samgene.'.quality.xls  ';#########added by dushuo
		my @clu_photos;
		push @clu_photos,"$sample_clu_dir/$samgenepdf.length.svg","$sample_clu_dir/$samgenepdf.length.png","$sample_clu_dir/$samgenepdf.length.txt";
		foreach my $clu_result(@clu_photos){
			my $clu_filename = basename($clu_result);
			my $clu_target = $up_clu_dir.'/'.$clu_filename;
			my $clu_shell = &_Upload($clu_result,$clu_target);
			print UPLOADSH $clu_shell;
		}
		undef @clu_photos;
		print UPLOADSH <<TGICLUP;
perl $config{"get_chosen_fa"} -fa $sample_clu_dir/$samgene -output $up_clu_dir/$samgene  -size -gap && rm $up_clu_dir/$samgene.list  && echo finish  upload  $sam-Unigen.fa at time $time

TGICLUP
		print CLEARSH <<TGICLRM;
rm  -r  $sample_clu_dir/asm_*   $sample_clu_dir/for_$sam.cluster.fa.* && echo rm  tgicl mid files at time $time 
rm $sample_clu_dir/cluster_and_single.fa         $sample_clu_dir/single.fa  $sample_clu_dir/$sam-Unigene.fa.quality.xls   && echo finish rm  $sam  cluster  mid files at time  $time 

TGICLRM
print VERSION 'perl  '.$config{svn_version}.'      '.$config{get_chosen_fa}.'  '.$config{replace_geneid}.'  '.$config{fa_quality}.' '.$config{barplot}.' '.$config{qsub_sge}.'  > '.$version_dir.'/sam_cluster.txt '."\n";
	}else{
		print STDERR 'warnings : no cluster of '.$sam."\n".'you can add tgicl_option=xxxx in lib file to  cluster or give the  clustered  fa file in lib file by add unigene= xxxx'."\n";
	}
}
my $quality_dir=$outdir.'/quality';
system "mkdir -p  $quality_dir" unless(-d $quality_dir);	
my $quality_sam_sh=$shdir.'/qua_sam.sh';
open SAMQUA,'>',$quality_sam_sh or die "can't open the sh of sam as quality";
#print SAMQUA 'perl '.$config{QCStatistics}.' '.$outdir.' '.$lib."\n";
if ($cat_qua_contig ne '') {
	print SAMQUA 'perl  '.$config{denovo_quality}.' -outdir '.$quality_dir.' -min 100  -type contig '.$cat_qua_contig."\n";
	my $contig_qua_shell = &_Upload("$quality_dir/contig_quality.xls","$up_as_fi_dir/1.Contig/Contig.quality.xls");
	print UPLOADSH $contig_qua_shell;
}
#######################################################added by dushuo to replace the three lines behind
if ($for_qua_file ne '') {
	print SAMQUA 'perl  '.$config{denovo_quality}.' -outdir '.$quality_dir.'  -type cluster '.$for_qua_file."\n";
}
close SAMQUA;

##########################################################################################
#############  聚类
##########################################################################################
my $clustersh=$shdir.'/cluster.sh';
my $clu_dir=$outdir.'/cluster';
system "mkdir -p  $clu_dir" unless(-d $clu_dir);
open CLUSH,'>',$clustersh or die "can't open the sh of cluster";
print CLUSH 'export LD_LIBRARY_PATH="'.$rna_lib.':$LD_LIBRARY_PATH"';
print CLUSH "\n";
my $gene='';
my $genepdf;
if(exists $option{$cluster}{unigene}){
	system "mkdir -p  $up_clu_dir " unless(-d $up_clu_dir);
	$option{$cluster}{unigene}=`perl $config{absolute_dir} $option{$cluster}{unigene}`;
	unless(-s $option{$cluster}{unigene}) {
		die 'warning  : '.$option{$cluster}{unigene}.' is not exists or it size is zero '."\n";
	}
	print CLUSH 'ln -s '.$option{$cluster}{unigene}.'  '.$clu_dir.'  && echo finish ln '.$clu_dir.'  at time '.$time."\n";
	$gene=basename($option{$cluster}{unigene});
	$genepdf=$gene;
	$genepdf=~s/\.fa$//;
	my @clusters;
	push @clusters,"$clu_dir/$gene","$clu_dir/$genepdf.length.svg","$clu_dir/$genepdf.length.png","$clu_dir/$genepdf.length.txt";
	foreach my $clu_result(@clusters){
		my $clu_filename = basename($clu_result);
		my $clu_target = $up_clu_dir.'/'.$clu_filename;
		my $clu_shell = &_Upload($clu_result,$clu_target);
		print UPLOADSH $clu_shell;
	}
	undef @clusters;
	$for_clu_file.=$clu_dir.'/'.$gene.'  ';
	$for_qua_file.=$clu_dir.'/'.$gene.'.quality.xls  ';###########added by dushuo###########
	print CLUSH 'perl '.$config{fa_quality}.' -len -Head -gap -N  -gc '.$clu_dir.'/'.$gene.' && echo finish statistics  All-Unigene.fa for draw pdf at time '.$time."\n";
	print CLUSH ' perl '.$config{barplot}.'   '.$clu_dir.'/'.$gene.'.quality.xls  '.$genepdf.'  -gap   && echo finish draw All-Unigene.pdf at time '.$time."\n"; 
	print CLUSH "\n";
	$unigene_file = $clu_dir.'/'.$gene;
}elsif (scalar(@samplename)==1){
	system "mkdir -p  $up_clu_dir " unless(-d $up_clu_dir);
	print CLUSH 'ln -s '.$outdir.'/assembly_fill/'.$samplename[0].'/cluster/*  '.$clu_dir.' && echo finish ln unigene at time '.$time."\n";
	$gene=$samplename[0].'-Unigene.fa';
	$genepdf=$gene;
	$genepdf=~s/\.fa$//;
	$unigene_file = $clu_dir.'/'.$gene;
}elsif (exists $option{$cluster}{tgicl_option}){
	system "mkdir -p  $up_clu_dir " unless(-d $up_clu_dir);
	$gene='All-Unigene.fa';
	$genepdf=$gene;
	$genepdf=~s/\.fa$//;
	$option{$cluster}{"unigene_len"} ||=200;
	print CLUSH <<TGICL;
cd $clu_dir && echo cd $clu_dir
cat  $for_clu_file  >   for_cluster.fa  &&  echo finish   get for_cluster.fa at time $time 
perl $config{tgicl} for_cluster.fa $option{$cluster}{tgicl_option}  && echo finish tgicl at $time 
cat asm_*/align > align   && echo get align at $time  
perl $config{"phrap.id.list.pl"} align align && echo get cluster id  at $time   
perl $config{"get_single.pl"}  align.cluster for_cluster.fa single  && echo get single  at $time 
cat asm_*/contigs > asm_cluster.fa  && echo get cluster sequences at $time
cat single.prefect.fa >> single.fa  && echo get singleton sequences at $time
cat asm_*/contigs single.fa > tgicl_cluster_and_single.fa && echo get cluster out file  at $time
cat align.cluster single.list > tgicl_cluster_and_single.fa.list && echo get tgicled unigenes list file at $time
/opt/blc/genome/bin/formatdb -p F -i tgicl_cluster_and_single.fa && echo formatdb tgicl_cluster_and_single.fa at $time
/opt/blc/genome/bin/blastall -p blastn -m 8 -e 1e-10 -F F -a 10 -d $clu_dir/tgicl_cluster_and_single.fa -i $clu_dir/tgicl_cluster_and_single.fa -o $clu_dir/all_vs_all.blast.m8 && echo blastn all to all at $time
perl $config{"cluster_for_coverage.pl"} $clu_dir/all_vs_all.blast.m8 $clu_dir && echo cluster based on coverage at $time
perl $config{"clusterOrSingle.pl"} $clu_dir All $option{$cluster}{"unigene_len"} && echo rename cluster and chose multiple cluster at $time
perl $config{"get_chosen_fa"} -fa All.all.fa -output All-Unigene.fa -len $option{$cluster}{"unigene_len"} && echo get $option{$cluster}{"unigene_len"} cluster fa  at $time
perl $config{"get_chosen_fa"} -fa cluster.all.fa -output All-Cluster.fa -len $option{$cluster}{"unigene_len"} && echo get $option{$cluster}{"unigene_len"} cluster fa at $time
perl $config{"get_chosen_fa"} -fa singleton.all.fa -output All-Single.fa -len $option{$cluster}{"unigene_len"} && echo get $option{$cluster}{"unigene_len"} singleton fa at $time
mkdir -p id_tmp
ln -s $clu_dir/cluster_and_single.fa.id.list id_tmp/cluster_and_single.fa.id.list
perl $config{"fa_quality"} -len -Head -gap -N  -gc   $gene && echo fa quality   at $time 
perl $config{"barplot"}   $gene.quality.xls   $genepdf  -gap  && echo plot   at $time 

TGICL
	$unigene_file = $clu_dir.'/'.$gene;
	$for_qua_file.=$clu_dir.'/'.$gene.'.quality.xls  ';###########added by dushuo###########

		print UPLOADSH <<TGICLUP;
perl $config{"get_chosen_fa"} -fa $clu_dir/$gene -output $up_clu_dir/$gene -size -gap && rm $up_clu_dir/$gene.list  && echo finish  upload  $gene at time $time

TGICLUP
		my @tgicls;
		push @tgicls,"$clu_dir/$genepdf.length.svg","$clu_dir/$genepdf.length.png","$clu_dir/$genepdf.length.txt";
		foreach my $tgi_result(@tgicls){
			my $tgi_filename = basename($tgi_result);
			my $tgi_target = $up_clu_dir.'/'.$tgi_filename;
			my $tgi_shell = &_Upload($tgi_result,$tgi_target);
			print UPLOADSH $tgi_shell;
		}
		undef @tgicls;
		my $tgi_list_shell = &_Upload("$clu_dir/cluster_and_single.fa.id.list","$up_clu_dir/$genepdf.id.xls");
		print UPLOADSH $tgi_list_shell;
		print CLEARSH <<TGICLRM;
rm  -r  $clu_dir/asm_*   $clu_dir/for_cluster.fa.* && echo rm  tgicl mid files at time $time 
rm $clu_dir/cluster_and_single.fa         $clu_dir/single.fa  $clu_dir/$gene.quality.xls   && echo finish rm $gene  cluster  mid files at time  $time 

TGICLRM

}else{
	print STDERR  'warnings :  no  Unigene  fa of '.$cluster."\n".'you can add tgicl_option=xxxx in lib file to  cluster or give the  clustered  fa file in lib file by add unigene= xxxx'."\n";
}
close CLUSH;
####################################################### quality	
my $quslity_sh=$shdir.'/quality.sh';
open QUA,'>',$quslity_sh or die "can't open the quality  sh ";
if ($for_qua_file ne '') {#############changed by dushuo
	print QUA 'perl  '.$config{denovo_quality}.' -outdir '.$quality_dir.'  -type Unigene '.$for_qua_file."\n";############changed by dushuo
	print UPLOADSH 'cp  '.$quality_dir.'/Unigene_quality.xls   '.$up_as_fi_dir.'/2.Unigene/Unigene.quality.xls  && echo finish upload Unigene.quality.xls at time '.$time."\n"; 
	my $uni_qua_shell = &_Upload("$quality_dir/Unigene_quality.xls","$up_as_fi_dir/2.Unigene/Unigene.quality.xls");
	print UPLOADSH $uni_qua_shell;
	print QUA 'perl '.$config{assemblyQuality}.' '.$quality_dir."\n";
	print UPLOADSH 'cp '.$quality_dir.'/assembly_statistic.xls '.$up_as_fi_dir.'/assembly_statistic.xls  && echo finish upload assembly_statistic.xls at time '.$time."\n";
	my $ass_sta_shell = &_Upload("$quality_dir/assembly_statistic.xls","$up_as_fi_dir/assembly_statistic.xls");
}
if ($cat_qua_fq ne '') {
#	print QUA 'perl '.$config{fq_num_bp}.' -out '.$quality_dir.'/Sequencing_output.xls  '.$cat_qua_fq."\n\n\n";
	#print UPLOADSH 'cp  '.$quality_dir.'/Sequencing_output.xls   '.$up_as_fi_dir.'  && echo finish upload Sequencing_output.xls at time '.$time."\n"; 
	my $bp_num_shell = &_Upload("$quality_dir/Sequencing_output.xls","$up_as_fi_dir/Sequencing_output.xls");
	print UPLOADSH $bp_num_shell;
}
my $who=`whoami`;
chomp $who;
my $report_time=`date +%y-%m-%d.%H`;
chomp $report_time;
$option{$cluster}{"denovo_info"} ||="/ifs2/BC_GAG/Project/RNA_denovo/xutong/denovo_info/$option{$cluster}{pro_name}_denovo_info.xls";
my $all_denovo_info_dir=dirname($option{$cluster}{"denovo_info"});
unless (-d  $all_denovo_info_dir) {
	if (-d "/ifs2/BC_GAG/") {
		$option{$cluster}{"denovo_info"}="/ifs2/BC_GAG/Project/RNA_denovo/xutong/denovo_info/$option{$cluster}{pro_name}_denovo_info.xls";
	}else{
		$option{$cluster}{"denovo_info"}="/ifshk1/BC_gag/PROJECT/RNA_denovo/xutong/data/$option{$cluster}{pro_name}_denovo_info.xls";
	}
}
my $all_denovo_info=$option{$cluster}{"denovo_info"};
my $ownerinfo="$who"."_$report_time";
my $table_head;
if(-f $all_denovo_info){
	$table_head = "n";
}else{
	$table_head = "y";
	system "mkdir -p $all_denovo_info && rmdir $all_denovo_info";
}
print QUA "perl $config{get_denovo_info} -lib $lib -fq $outdir/fq -dir $outdir -outfile $all_denovo_info -table_head $table_head -ownerinfo $ownerinfo  2> /dev/null\n";
print QUA  'echo  finish  quality at `date +%y-%m-%d.%H:%M:%S` '."\n\n\n";
close QUA;

##########################################################################################
#############  注释
##########################################################################################
my $gene_name=$gene;
$gene_name=~s/\.fa$// if ($gene ne '');
my $blast_xls_fa_name;
my $blast_sub_out_name;
my $up_blast_dir=$upload_dir.'/annotation';
my $annotdir=$outdir.'/annot';
my $godir=$annotdir.'/go';
my $go_sh=$shdir.'/go.sh';
my $cds_sh=$shdir.'/cds.sh';
my $fs_sign_sh=$shdir.'/sign.sh';
########################################annot.pl --> by dushuo
my $annot_pl="$shdir/annot.pl";
open ANNO,'>',$annot_pl;
print ANNO<<LYSIS;
#!/usr/bin/perl
use strict;
use lib "/opt/rocks/lib/perl5/5.10.1";
use Thread 'async';

my \$blast_go = async {
LYSIS
if ((exists $option{$cluster}{blast_option}) && ($gene ne '')) {

	print ANNO<<LYSIS;
	my \$step1=1;
	\$step1=system ("perl $shdir/blast_annot.pl");
LYSIS

	system "mkdir -p  $annotdir" unless (-d $annotdir);
	my $databaseannotdir=$annotdir.'/database';
	system "mkdir -p  $databaseannotdir" unless (-d $databaseannotdir);
	my $blast_xls_outdir=$databaseannotdir.'/'.$gene_name.'/gene-annotation';
	my $blast_xls_head=$blast_xls_outdir.'/blast_tab_head  ';
	my $blast_xls_head_m8 = $blast_xls_outdir.'/blast_tab_head.m8';
	$blast_xls_fa_name=$blast_xls_outdir.'/'.$gene;
	$blast_sub_out_name=$blast_xls_outdir.'/'.$gene_name.'.fa.cut/'.$gene_name.'.fa';
	system "mkdir -p  $up_blast_dir" unless(-d $up_blast_dir);
	my $blast_xls='';
############### --> dushuo
	my $blast_annot=$shdir."/blast_annot.pl";
	open ANNOT,'>',$blast_annot;
##############
print ANNOT<<BANNOT;
#!/usr/bin/perl
use strict;
use lib "/opt/rocks/lib/perl5/5.10.1";
use Thread 'async';

BANNOT
	##############  blast#####changed by dushuo 2010-11-24
	my $blast_cwd.='system ("perl '.$config{search_database}.'  '.$option{$cluster}{blast_option}.'  -input '.$clu_dir.'/'.$gene.'   -outdir  '.$databaseannotdir.'  -path1 '.$path1.' -path2 '.$path2." -shdir $shdir ".' && echo finish blast at time '.$time."\");\n"; 
	$blast_cwd.='system ("cp '.$databaseannotdir.'/blast_db_info.txt  '.$version_dir.'/blast_db_info.txt ");'."\n\n";
	$cwdmv='system ("find  '.$shdir.' -maxdepth 1 -name \"'.$gene.'.blast.*sh.*\"  -exec mv {} '.$oldshdir.'\";\"");';
	print ANNOT $cwdmv."\n";
	print ANNOT $blast_cwd."\n";
	my $judge='';
############################### V structure changed by dushuo 2010-11-24
	if ($option{$cluster}{blast_option}=~/-cog/) {
		print ANNOT<<COG;
my \$cog = async {
	my \$step1=1;
	\$step1=system ("perl $config{qsub_sge}  --reqsub  --getmem  --queue $queue1 --subprjctid $option{$cluster}{sub_prjct_id} --resource vf=1G --jobprefix blast --lines 1 --interval 300 $shdir/$gene.blast.cog.sh ");
	my \$step2=1;
	\$step2=system("perl $config{finish}  -indir $shdir/$gene.blast.cog.sh.*.qsub") if (\$step1 eq 0);
	my \$step3=1;
	\$step3=system("sh $shdir/cat_blast.cog.$gene.sh") if (\$step2 eq 0);
	if(\$step3 eq 0){return 0}else{return 1}
};
COG
		$judge.="\$cog->join() == 0 && ";
	}
	if ($option{$cluster}{blast_option}=~/-kegg/) {
		print ANNOT<<KEGG;
my \$kegg = async {
	my \$step1=1;
	\$step1=system ("perl $config{qsub_sge}  --reqsub  --getmem  --queue $queue1 --subprjctid $option{$cluster}{sub_prjct_id} --resource vf=1G --jobprefix blast --lines 1 --interval 300 $shdir/$gene.blast.kegg.sh ");
	my \$step2=1;
	\$step2=system("perl $config{finish}  -indir $shdir/$gene.blast.kegg.sh.*.qsub") if (\$step1 eq 0);
	my \$step3=1;
	\$step3=system("sh $shdir/cat_blast.kegg.$gene.sh") if (\$step2 eq 0);
	if(\$step3 eq 0){return 0}else{return 1}
};
KEGG
		$judge.="\$kegg->join() == 0 && ";
	}
	if ($option{$cluster}{blast_option}=~/-swissprot/) {
		print ANNOT<<SWISS;
my \$swissprot = async {
	my \$step1=1;
	\$step1=system ("perl $config{qsub_sge}  --reqsub  --getmem  --queue $queue1 --subprjctid $option{$cluster}{sub_prjct_id} --resource vf=1G --jobprefix blast --lines 1 --interval 300 $shdir/$gene.blast.swissprot.sh ");
	my \$step2=1;
	\$step2=system("perl $config{finish}  -indir $shdir/$gene.blast.swissprot.sh.*.qsub") if (\$step1 eq 0);
	my \$step3=1;
	\$step3=system("sh $shdir/cat_blast.swiss.$gene.sh") if (\$step2 eq 0);
	if(\$step3 eq 0){return 0}else{return 1}
};
SWISS
		$judge.="\$swissprot->join() == 0 && ";
	}
############################# A
	if (($option{$cluster}{blast_option}=~/-nr/) || ($option{$cluster}{blast_option}=~/-userdb/) ){
		print ANNOT<<NR;
my \$nr = async {
	my \$step1=1;
	\$step1=system ("perl $config{qsub_sge}  --reqsub  --getmem  --queue $queue1 --subprjctid $option{$cluster}{sub_prjct_id} --resource vf=2.5G --jobprefix blast --lines 1 --interval 300 $shdir/$gene.blast.nr.sh ");
	my \$step2=1;
	\$step2=system("perl $config{finish}  -indir $shdir/$gene.blast.nr.sh.*.qsub") if (\$step1 eq 0);
	my \$step3=1;
	\$step3=system("sh $shdir/cat_blast.nr.$gene.sh") if (\$step2 eq 0);
	if(\$step3 eq 0){return 0}else{return 1}
};
NR
		$judge.="\$nr->join() == 0 && ";
	}
	if ($option{$cluster}{blast_option}=~/-nt/){
		print ANNOT<<NT;
my \$nt = async {
	my \$step1=1;
	\$step1=system ("perl $config{qsub_sge}  --reqsub  --getmem  --queue $queue1 --subprjctid $option{$cluster}{sub_prjct_id} --resource vf=6.5G --jobprefix blast --lines 1 --interval 300 $shdir/$gene.blast.nt.sh ");
	my \$step2=1;
	\$step2=system("perl $config{finish}  -indir $shdir/$gene.blast.nt.sh.*.qsub") if (\$step1 eq 0);
	my \$step3=1;
	\$step3=system("sh $shdir/cat_blast.nt.$gene.sh") if (\$step2 eq 0);
	if(\$step3 eq 0){return 0}else{return 1}
};
NT
		$judge.="\$nt->join() == 0 && ";
	}
	$judge=~s/&&\s*$//;
	print ANNOT<<BLAST;
if ($judge) {
	print "Finish : blast nr cog kegg swissprot nt!";
	exit 0;
}else{
	print "Warn : blast in one of nr cog kegg swissprot nt went wrong!";
	exit 1;
}
BLAST
	close ANNOT;

	print CLEARSH ' rm -r '.$blast_xls_fa_name.'.cut  && echo finish rm blast cut result  at time '.$time."\n";
	if ($option{$cluster}{blast_option}=~/-cog/) {
		system "mkdir -p  $up_blast_dir/COG " unless ( -d  "$up_blast_dir/COG" );
		$blast_xls.=$blast_xls_fa_name.'.blast.cog.xls   ';
		print UPLOADSH 'cat   '.$blast_xls_head_m8.'  '.$blast_xls_fa_name.'.blast.cog.xls    >  '.$up_blast_dir.'/COG/'.$gene.'.blast.cog.xls   && echo finish upload blast.cog.xls  at time '.$time."\n";	
		my @cogs;
		push @cogs,"$blast_xls_fa_name.cog.class.annot.xls","$blast_xls_fa_name.cog.pdf","$blast_xls_fa_name.cog.png","$blast_xls_fa_name.cog.gene.annot.xls";
		foreach my $cog_result(@cogs){
			my $cog_filename = basename($cog_result);
			my $target = $up_blast_dir.'/COG/'.$cog_filename;
			my $shell_script = &_Upload($cog_result,$target);
			print UPLOADSH $shell_script;
		}
		undef @cogs;
	}
	if ($option{$cluster}{blast_option}=~/-kegg/) {
		system "mkdir -p  $up_blast_dir/KEGG" unless(-d "$up_blast_dir/KEGG" );
		$blast_xls.=$blast_xls_fa_name.'.blast.kegg.xls   ';
		print UPLOADSH 'cat '.$blast_xls_head_m8.' '.$blast_xls_fa_name.'.blast.kegg.xls > '.$up_blast_dir.'/KEGG/'.$gene.'.blast.kegg.xls   && echo finish upload blast.kegg.xls  at time '.$time."\n";
		my @keggs;
		push @keggs,"$blast_xls_fa_name.htm","$blast_xls_fa_name.ko","$blast_xls_fa_name.path","$blast_xls_fa_name\_map";
		foreach my $kegg_result (@keggs){
			my $kegg_filename = basename($kegg_result);
			my $target = $up_blast_dir.'/KEGG/'.$kegg_filename;
			my $shell_script = &_Upload($kegg_result,$target);
			print UPLOADSH $shell_script;
		}
		undef @keggs;
		print CLEARSH 'rm '.$blast_xls_fa_name.'.blast.kegg  && echo finish rm blast kegg m0 result file at time '.$time."\n";
	}
	if ($option{$cluster}{blast_option}=~/-nr/) {
		system "mkdir -p  $up_blast_dir/Nr" unless(-d "$up_blast_dir/Nr" );
		$blast_xls.=$blast_xls_fa_name.'.blast.Nr.xls   ';
		print UPLOADSH 'cat   '.$blast_xls_head.'  '.$blast_xls_fa_name.'.blast.Nr.xls   >  '.$up_blast_dir.'/Nr/'.$gene.'.blast.Nr.xls  && echo finish upload blast.Nr.xls at time '.$time."\n";
		print UPLOADSH "cp $blast_xls_outdir/*_statistic.xls $up_blast_dir/Nr/
cp $blast_xls_outdir/class_statistics.png $up_blast_dir/Nr/ && echo finish upload nr statistics files at time $time\n"; 
		print CLEARSH 'rm '.$blast_xls_fa_name.'.blast.Nr* && echo finish rm blast nr m0 files at time '.$time."\n";
	}
	if ($option{$cluster}{blast_option}=~/-swissprot/) {
		system "mkdir -p  $up_blast_dir/Swissprot" unless(-d "$up_blast_dir/Swissprot" );
		$blast_xls.=$blast_xls_fa_name.'.blast.Swissprot.xls   ';
		print UPLOADSH 'cat   '.$blast_xls_head_m8.'  '.$blast_xls_fa_name.'.blast.Swissprot.xls   >  '.$up_blast_dir.'/Swissprot/'.$gene.'.blast.Swissprot.xls   &&  echo finish  upload blast.Swissprot.xls at time '.$time."\n"; 
		print CLEARSH 'rm '.$blast_xls_fa_name.'.blast.Swissprot && echo finish rm blast swissprot m0 files at time '.$time."\n";
	}
	if ($option{$cluster}{blast_option}=~/-userdb\s(\S+)/) {
		my $usdb=basename($1);
		system "mkdir -p  $up_blast_dir/$usdb" unless(-d "$up_blast_dir/$usdb" );
		$blast_xls.=$blast_xls_fa_name.'.blast.'.$usdb.'.xls   ';
		print UPLOADSH 'cat   '.$blast_xls_head.'  '.$blast_xls_fa_name.'.blast.'.$usdb.'.xls   >  '.$up_blast_dir.'/'.$usdb.'/'.$gene.'.blast.'.$usdb.'.xls   &&  echo finish  upload blast.'.$usdb.'.xls at time '.$time."\n"; 
		print CLEARSH 'rm '.$blast_xls_fa_name.'.blast.'.$usdb.' && echo finish rm blast '.$usdb.' m0 files at time '.$time."\n";
	}
	if ($option{$cluster}{blast_option}=~/-nt/) {
		system "mkdir -p  $up_blast_dir/Nt" unless(-d "$up_blast_dir/Nt" );
		$blast_xls.=$blast_xls_fa_name.'.blast.Nt.xls   ';
		print UPLOADSH 'cat   '.$blast_xls_head.'  '.$blast_xls_fa_name.'.blast.Nt.xls   >  '.$up_blast_dir.'/Nt/'.$gene.'.blast.Nt.xls  && echo finish upload blast.Nt.xls at time '.$time."\n"; 
		print CLEARSH 'rm '.$blast_xls_fa_name.'.blast.Nt && echo finish rm blast nt m0 files at time '.$time."\n";
	}

	############# go
	if (exists $option{$cluster}{nr2go_option}) {
		my $nr_file=$blast_xls_fa_name.'.blast.Nr';
		my $m7file=$godir.'/'.$gene.'.blast.Nr.xml';
		system "mkdir -p  $godir" unless (-d $godir);

		open GOSH,'>',$go_sh or die "can't open the sh file of go $go_sh";
		my $gooption=$1 if ($option{$cluster}{nr2go_option} =~/(-gvol\s*\d+)/);
		print CLEARSH 'rm '.$m7file.'  && echo finish rm nr2go m7 file at time '.$time."\n";
		print GOSH 'cat '.$blast_sub_out_name.'.*.blast.Nr.xml.annot > '.$godir.'/'.$gene.'.blast.Nr.xml.annot && echo finish cat annot files at time '.$time."\n";
		print GOSH 'perl '.$config{annot2wego}.'  -i  '.$m7file.'.annot -outdir '.$godir.'/    && echo finish get wego file at time '.$time."\n";
		print GOSH 'perl '.$config{drawGO}.' -gglist  '.$m7file.'.wego  -output  '.$godir.'/'.$gene.'.GO  -go  '.$config{go_class}.'    && echo finish draw svg at time '.$time."\n";
		#print GOSH 'perl '.$config{changsvgsize}.'  '.$godir.'/'.$gene.'.GO.svg  50  && echo  finish changsvgsize at time '.$time."\n";
		print GOSH $config{"java"}.' -Djava.awt.headless=true  -jar '.$config{"batik-rasterizer.jar"}.' -m image/png '.$godir.'/'.$gene.'.GO.svg && echo finish  svg  to png at time '.$time."\n";
		my @tmp = split(/\-Unigene\.fa/, $gene);
		print GOSH 'perl '.$config{"annot_statistic"};
		if($option{$cluster}{blast_option}=~/-nr/){
			print GOSH ' -nr '.$blast_xls_fa_name.'.blast.Nr.xls';
		}
		if($option{$cluster}{blast_option}=~/-nt/){
			print GOSH ' -nt '.$blast_xls_fa_name.'.blast.Nt.xls';
		}
		if($option{$cluster}{blast_option}=~/-swissprot/){
			print GOSH ' -swissprot '.$blast_xls_fa_name.'.blast.Swissprot.xls';
		}
		if($option{$cluster}{blast_option}=~/-kegg/){
			print GOSH ' -kegg '.$blast_xls_fa_name.'.blast.kegg.xls';
		}
		if($option{$cluster}{blast_option}=~/-cog/){
			print GOSH ' -cog '.$blast_xls_fa_name.'.blast.cog.xls';
		}
		print GOSH ' -go '.$m7file.'.wego -sample '.$gene.' -out '.$outdir.'/quality/annotation_statistic.xls'."\n";
		#print GOSH ' perl '.$config{svg2xxx}.'   '.$godir.'/'.$gene.'.GO.svg     && echo finish  svg  to png at time '.$time."\n";
		close GOSH;

		my $godb=<<GODB;
blast2go = $config{blast2go}
b2gPipe_properties = $config{b2gPipe_properties}
go_class = $config{go_class}
go_alias = $config{go_alias}
GODB

		open GODBFILE,'>',$version_dir.'/go_db.txt' or die "can't open the go db file info";
		print GODBFILE $godb;
		close GODBFILE;

		open WEGOHEAD,'>',$godir.'/wego.head' or die "can't open the head of wego ";
		print WEGOHEAD 'geneID'."\t".'GO'."\n";
		close WEGOHEAD;

		system "mkdir -p  $up_blast_dir/GO " unless (-d "$up_blast_dir/GO" );
		print UPLOADSH 'cat   '.$godir.'/wego.head '.$m7file.'.wego  > '.$up_blast_dir.'/GO/'.$gene.'.gene2GO.xls   && echo finish upload gene2GO.xls at time  '.$time."\n";
		print UPLOADSH 'cp '.$outdir.'/quality/annotation_statistic.xls '.$up_blast_dir.'/annotation_statistic.xls  && echo finish upload annotation_statistic.xls at time  '.$time."\n";
		my $shell_scri = &_Upload("$godir/$gene.GO.xls","$up_blast_dir/GO/$gene.GO2gene.xls");
		print UPLOADSH $shell_scri;
		my @gos;
		push @gos,"$godir/$gene.GO.svg","$godir/$gene.GO.png";
		foreach my $go_result(@gos){
			my $go_filename = basename($go_result);
			my $target = $up_blast_dir.'/GO/'.$go_filename;
			my $shell_script = &_Upload($go_result,$target);
			print UPLOADSH $shell_script;
		}
	}
	##########################################################################################
	############# nr ESTscan
	##########################################################################################
	my $cdsdir=$outdir.'/CDS';
	system "mkdir -p  $cdsdir " unless(-d $cdsdir);
	my $estdir=$cdsdir.'/ESTscan';
	system "mkdir -p  $estdir " unless(-d $estdir);

	$option{$cluster}{"cds_len"} ||=100;
	unless (exists $option{$cluster}{"nrEST_option"}) {
		$option{$cluster}{"nrEST_option"}='  -l '.$option{$cluster}{"cds_len"};
	}
	my $nrEST_option_l=$option{$cluster}{"nrEST_option"};

	my $cwd_cat_xls="";##by chuanjun
	if ($option{$cluster}{blast_option}=~/-userdb\s(\S+)/) {
		my $usdb=basename($1);
		$cwd_cat_xls.="$blast_xls_fa_name.blast.$usdb.xls,";
	}
	if ($option{$cluster}{blast_option}=~/-nr/) {
		$cwd_cat_xls.="$blast_xls_fa_name.blast.Nr.xls,";
	}
	if ($option{$cluster}{blast_option}=~/-swissprot/) {
		$cwd_cat_xls.="$blast_xls_fa_name.blast.Swissprot.xls,";
	}
	if ($option{$cluster}{blast_option}=~/-kegg/) {
		$cwd_cat_xls.="$blast_xls_fa_name.blast.kegg.xls,";
	}
	if ($option{$cluster}{blast_option}=~/-cog/) {
		$cwd_cat_xls.="$blast_xls_fa_name.blast.cog.xls";
	}
	open CDSSH,'>',$cds_sh or die "can't open the sh of get cds form blast  and ESTscan $cds_sh";
	print CDSSH <<TEMPSH;
perl $config{"get_cds_blast.pl"} -fa  $clu_dir/$gene -xls $cwd_cat_xls  -out $cdsdir/$gene_name.blast  -L $option{$cluster}{"cds_len"} 
rm -r $estdir/* 
mv $cdsdir/$gene_name.blast.mrna.fa  $estdir/mrna.seq 
cp $cdsdir/$gene_name.conf  $estdir 
perl $config{"prepare_data"}   -e $cdsdir/$gene_name.conf 
perl $config{"build_model"}     $cdsdir/$gene_name.conf 
$config{estscan} $cdsdir/$gene_name.blast.no.fa -o $cdsdir/$gene_name.ESTscan.cds.fa.score -t $cdsdir/$gene_name.ESTscan.protein.fa.score   -M  $estdir/Matrices/*.smat 
perl $config{"clear.score.pl"}  $cdsdir/$gene_name.ESTscan.cds.fa.score $cdsdir/$gene_name.ESTscan.cds.fa -debug  
perl $config{"clear.score.pl"}  $cdsdir/$gene_name.ESTscan.protein.fa.score  $cdsdir/$gene_name.ESTscan.protein.fa -debug 
perl $config{"fa_quality"}  -len -Head -gap -N  -gc  $cdsdir/$gene_name.blast.cds.fa   $cdsdir/$gene_name.ESTscan.cds.fa 
perl $config{"fa_quality"}  -len -Head  $cdsdir/$gene_name.blast.protein.fa  $cdsdir/$gene_name.ESTscan.protein.fa 
perl $config{"barplot"}  $cdsdir/$gene_name.blast.cds.fa.quality.xls  $gene_name.blast.cds.fa -gap 
perl $config{"barplot"}  $cdsdir/$gene_name.ESTscan.cds.fa.quality.xls  $gene_name.ESTscan.cds.fa -gap 
perl $config{"barplot"}  $cdsdir/$gene_name.blast.protein.fa.quality.xls  $gene_name.blast.protein.fa 
perl $config{"barplot"}  $cdsdir/$gene_name.ESTscan.protein.fa.quality.xls  $gene_name.ESTscan.protein.fa 
TEMPSH
	close CDSSH;
	open ESTCONF,'>',$cdsdir.'/'.$gene_name.'.conf' or die "can't open the conf of ESTscan $cdsdir/$gene_name.conf";
	print ESTCONF <<TEMPCONG;

################################################################################
#
# Parameters for $gene_name
# (use PERL syntax!)
#
\$organism = \"$gene_name\";\
\$hightaxo = \"\";
\$dbfiles =\"\";
\$ugdata = \"\";
\$estdata = \"\";
\$datadir = \"$estdir\";
\$nb_isochores = 2;
\$tuplesize = 6;
\$minmask = 30;
#
# End of File
#
################################################################################

TEMPCONG
	close ESTCONF;

	print CLEARSH 'rm -r '.$estdir.'/*  && echo finish rm ESTsan cds mid files at time '.$time."\n";

	my $up_cds_dir=$upload_dir.'/CDS';
	system "mkdir -p  $up_cds_dir" unless(-d $up_cds_dir);
	my @cdss;
	push @cdss,"$cdsdir/$gene_name.blast.cds.fa","$cdsdir/$gene_name.blast.protein.fa","$cdsdir/$gene_name.ESTscan.cds.fa","$cdsdir/$gene_name.ESTscan.protein.fa","$cdsdir/$gene_name.blast.cds.fa.length.svg","$cdsdir/$gene_name.blast.cds.fa.length.png","$cdsdir/$gene_name.blast.cds.fa.length.txt","$cdsdir/$gene_name.blast.protein.fa.length.svg","$cdsdir/$gene_name.blast.protein.fa.length.png","$cdsdir/$gene_name.blast.protein.fa.length.txt","$cdsdir/$gene_name.ESTscan.cds.fa.length.svg","$cdsdir/$gene_name.ESTscan.cds.fa.length.png","$cdsdir/$gene_name.ESTscan.cds.fa.length.txt","$cdsdir/$gene_name.ESTscan.protein.fa.length.svg","$cdsdir/$gene_name.ESTscan.protein.fa.length.png","$cdsdir/$gene_name.ESTscan.protein.fa.length.txt";
	foreach my $cds_result(@cdss){
		my $cds_filename = basename($cds_result);
		my $target = $up_cds_dir.'/'.$cds_filename;
		my $shell_script = &_Upload($cds_result,$target);
		print UPLOADSH $shell_script;
	}
	print UPLOADSH 'cp ',$config{"win_blast"},'/* ',$up_cds_dir ,' &&echo finish cp win_blast at time ',$time,"\n";
	##########################################################################################
	############# gene orientation
	##########################################################################################
	open SIGNFASH,'>',$fs_sign_sh or die "can't open the sh of orientation $fs_sign_sh";
	print SIGNFASH 'perl '.$config{sign};
	print SIGNFASH '  -ESTscan  '.$cdsdir.'/'.$gene_name.'.ESTscan.cds.fa  -fa  '.$clu_dir.'/'.$gene.'  -outdir '.$clu_dir;
	if ($option{$cluster}{blast_option}=~/-userdb\s(\S+)/) {
		my $usdb=basename($1);
		print SIGNFASH ' -userdb  '.$blast_xls_fa_name.'.blast.'.$usdb.'.xls  ';
	}
	if ($option{$cluster}{blast_option}=~/-nr/) {
		print SIGNFASH ' -nr '.$blast_xls_fa_name.'.blast.Nr.xls  ';
	}
	if ($option{$cluster}{blast_option}=~/-swissprot/) {
		print SIGNFASH ' -swissprot '.$blast_xls_fa_name.'.blast.Swissprot.xls  ';
	}
	if ($option{$cluster}{blast_option}=~/-cog/) {
		print SIGNFASH ' -cog '.$blast_xls_fa_name.'.blast.cog.xls   ';
	}
	if ($option{$cluster}{blast_option}=~/-kegg/) {
		print SIGNFASH ' -kegg '.$blast_xls_fa_name.'.blast.kegg.xls   ';
	}	
	print SIGNFASH '  && echo finish get the sign  from blast and ESTscan at time '.$time."\n";
	print SIGNFASH "\n";
	close SIGNFASH;
	my @signs;
	push @signs,"$clu_dir/$gene_name.orientation.xls","$clu_dir/$gene_name.5-3.fa","$clu_dir/$gene_name.no_orientation.fa";
	foreach my $sign_result(@signs){
		my $sign_filename = basename($sign_result);
		my $target = $up_clu_dir.'/'.$sign_filename;
		my $shell_scri = &_Upload($sign_result,$target);
		print UPLOADSH $shell_scri;
	}
}

##########################################################################################
####### 差异分析
##########################################################################################
my $bwt_sh=$shdir.'/bwt.sh';
my $soap_sh=$shdir.'/soap.sh';
my $exp_flag;
if(exists $option{$cluster}{exp_option}){
	if($option{$cluster}{exp_option} eq "rpkm"){
		$exp_flag = "rpkm";
	}
	elsif($option{$cluster}{exp_option} eq "fpkm"){
		$exp_flag = "fpkm";
	}
	else{
		die "No such algorithm as $option{$cluster}{exp_option}. Please check it!\n";
	}
}
my $exp_sh="$shdir/$exp_flag.sh";
my $pvalue_sh=$shdir.'/pvalue.sh';
my $functional_sh=$shdir.'/functional.sh';
#my $snp_sh=$shdir.'/first_snp.sh';
my $coverage_sh=$shdir.'/coverage.sh';

if (exists $option{$cluster}{soap}) {
	my $soap_dir=$outdir.'/soap';
	system "mkdir -p  $soap_dir" unless (-d $soap_dir);
	  #############  2bwt
#		my $bwt=$option{$cluster}{soap};
#		$bwt.='bwt';
		open BWT,'>',$bwt_sh or die "can't open the 2bwt sh $bwt_sh";
		if($gene ne '' ){
			print BWT ' ln -s '.$clu_dir.'/'.$gene.'  '.$soap_dir.'/'.$gene.' && echo finish get fa for soap at time '.$time."\n";
			print BWT $config{bwt}.'    '.$soap_dir.'/'.$gene.'  && echo finish bwt at time '.$time."\n";
			print BWT 'perl '.$config{fa_quality}.'    '.$soap_dir.'/'.$gene.' -len       && echo finish get len of gene at time '.$time."\n";
		}else{
			print STDERR 'warnings : no fa to bwt for soap '."\n".'you can give the  clustered  fa file in lib file by add unigene= xxxx'."\n";
		}
		close BWT;
		print CLEARSH 'rm '.$soap_dir.'/'.$gene.'.index.* && echo finish rm bwt index file'."\n";
	  #############  soap
		open SOAPSH,'>',$soap_sh or die "can't open the soap sh $soap_sh";
		my $soaplist=$soap_dir.'/soapfile.list';
		open SOAPLIST,'>',$soaplist or die "can't open the soap file list $soaplist ";
		my $sam_cov="$shdir/CoverSH";
		system "mkdir -p  $sam_cov" unless (-d $sam_cov);
		open TOTALCOVER,'>',$coverage_sh or die "Can't open the sh of soap coverage ";
		system "mkdir -p  $soap_dir/cover " unless (-d "$soap_dir/cover" );
		my $soap_out_suffix1=".PESoap";
		my $soap_out_suffix2=".PESoapSingle";
		if ($option{$cluster}{soap}=~/gz/) {
			$soap_out_suffix1.=".gz";
			$soap_out_suffix2.=".gz";
		}
		foreach my $sam_soap (@samplename) {
			my @fq=split /\n/,$option{$sam_soap}{fq};
			my $fq_num=@fq;
			my $i=0;
			my $sam_soap_num=1;
			print TOTALCOVER "sh $sam_cov/$sam_soap.cover.sh\n";
			open COVERSH,'>',"$sam_cov/$sam_soap.cover.sh" or die "Can't open the sh of soap coverage ";
			open SAMSOAPLIST,'>',$soap_dir.'/'.$sam_soap.'.soap.list' or die "can't open the soap list of $sam_soap ";
			while ($i<@fq) {
				my $soap=$option{$cluster}{soap};
				print SOAPSH $config{$soap};
				my $tag_pe=1;
				if ($fq[$i]=~/^q1\s*=\s*(\S+)/) {
					print SOAPSH '  -a '.$1;
					$i++;
					if ($fq[$i]=~/^q2\s*=\s*(\S+)/) {
						print SOAPSH '  -b '.$1;
					}else{
						print STDERR 'warnings :   error in lib of '.$sam_soap."\n";
					}			
				}elsif($fq[$i]=~/^q\s*=\s*(\S+)/){
					print SOAPSH '  -a '.$1;
					$tag_pe=0;
				}else{
					print STDERR 'warnings : error in lib of '.$sam_soap."\n";
				}
				$i++;
				print SOAPSH ' -D '.$soap_dir.'/'.$gene.'.index ';
				unless (exists $option{$cluster}{soap_option}) {
					$option{$cluster}{soap_option}=' -m 0 -x 1000 -s 40 -l 35 -v 3 ';
				}
				print SOAPSH $option{$cluster}{soap_option};
				if ($sam_soap_num>1) {
					if ($tag_pe==1) {
						print SOAPSH '  -o '.$soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix1.'  -2  '.$soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix2.'  ';
						print SOAPLIST $soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix1.' '.$soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix2."\t";
						print SAMSOAPLIST $soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix1."\n";
						print SAMSOAPLIST $soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix2."\n";
					}else{
						print SOAPSH '  -o '.$soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix2.'  ';
						print SOAPLIST $soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix2."\t";
						print SAMSOAPLIST $soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix2."\n";
					}
					print CLEARSH 'rm '.$soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix1.'  '.$soap_dir.'/'.$sam_soap.'_'.$sam_soap_num.$soap_out_suffix2.'  && echo finish rm '.$sam_soap.' soap '."\n";
				}else{
					if ($tag_pe==1) {
						print SOAPSH '  -o '.$soap_dir.'/'.$sam_soap.$soap_out_suffix1.'  -2  '.$soap_dir.'/'.$sam_soap.$soap_out_suffix2.' ';
						print SOAPLIST $soap_dir.'/'.$sam_soap.$soap_out_suffix1.' '.$soap_dir.'/'.$sam_soap.$soap_out_suffix2."\t";
						print SAMSOAPLIST $soap_dir.'/'.$sam_soap.$soap_out_suffix1."\n";
						print SAMSOAPLIST $soap_dir.'/'.$sam_soap.$soap_out_suffix2."\n";
					}else{
						print SOAPSH '  -o    '.$soap_dir.'/'.$sam_soap.$soap_out_suffix2.' ';
						print SOAPLIST $soap_dir.'/'.$sam_soap.$soap_out_suffix2."\t";
						print SAMSOAPLIST $soap_dir.'/'.$sam_soap.$soap_out_suffix2."\n";
					}
					print CLEARSH 'rm '.$soap_dir.'/'.$sam_soap.$soap_out_suffix1.' '.$soap_dir.'/'.$sam_soap.$soap_out_suffix2.'  && echo finish rm '.$sam_soap.' soap '."\n";
				}
				print SOAPSH  ' && echo finish  soap at time '.$time."\n";
				$sam_soap_num++;
			}	
			print SOAPLIST "\n";
			close SAMSOAPLIST;
		 ############# coverage
			print COVERSH 'perl '.$config{ReadsRandomInGene}.' '.$soap_dir.'/'.$gene.'.quality.xls '.$soap_dir.'/'.$sam_soap.'.soap.list '.$soap_dir.'/cover/'.$gene_name.'.'.$sam_soap.' && echo finish ReadsRandomInGene  at time '.$time."\n";
			print COVERSH $config{"java"}.' -Djava.awt.headless=true -jar '.$config{"batik-rasterizer.jar"}.' '.$soap_dir.'/cover/'.$gene_name.'.'.$sam_soap.'.ReadsRandom.svg && echo svg to png  at time '.$time."\n";
			print COVERSH 'perl '.$config{soapCoverage4RNAdenovo}.' -soap '.$soap_dir.'/'.$sam_soap.'.soap.list -unigene '.$unigene_file.' -outfile '.$soap_dir.'/cover/'.$gene_name.'.'.$sam_soap.'.Coverage.xls -verbose && echo finish Soap Coverage at time '.$time."\n";
			my @covers;
			push @covers,"$soap_dir/cover/$gene_name.$sam_soap.ReadsRandom.svg","$soap_dir/cover/$gene_name.$sam_soap.ReadsRandom.png","$soap_dir/cover/$gene_name.$sam_soap.Coverage.xls";
			foreach my $cover_result(@covers){
				my $cover_filename = basename($cover_result);
				my $cover_target = $up_clu_dir.'/'.$cover_filename;
				my $cover_shell = &_Upload($cover_result,$cover_target);
				print UPLOADSH $cover_shell;
			}
			close COVERSH;
		}
		close SOAPSH;
		close SOAPLIST;

	  ############# exp
		my $diff_gene=$soap_dir.'/'.$gene;
		my $genelen=$soap_dir.'/'.$gene.'.quality.xls';
		my $diff_dir=$outdir.'/diff';
		system "mkdir -p  $diff_dir" unless (-d $diff_dir);
		open EXP,'>',$exp_sh or die "can't open the $exp_flag sh $exp_sh";
		print EXP 'perl  '.$config{$exp_flag}.'  -len   '.$genelen.'  -list  '.$soaplist.'  -out '.$diff_dir.'/'.$exp_flag.'.xls  && echo finish get '.$exp_flag.' at time '.$time."\n";
		print EXP 'perl '.$config{get_sam_by_exp_new}.' -exp '.$diff_dir.'/'.$exp_flag.'.xls   -output '.$diff_dir.'/'.$exp_flag.'_annot.xls   ';
		if (exists $option{$cluster}{blast_option}) {
			if ($option{$cluster}{blast_option}=~/-userdb\s(\S+)/) {
				my $usdb=basename($1);
				print EXP ' -userdb  '.$blast_xls_fa_name.'.blast.'.$usdb.'.xls  ';
			}
			if ($option{$cluster}{blast_option}=~/-nr/) {
				print EXP ' -nr '.$blast_xls_fa_name.'.blast.Nr.xls  ';
			}
			if ($option{$cluster}{blast_option}=~/-swissprot/) {
				print EXP ' -swissprot '.$blast_xls_fa_name.'.blast.Swissprot.xls  ';
			}
			if ($option{$cluster}{blast_option}=~/-cog/) {
				print EXP ' -cog '.$blast_xls_fa_name.'.cog.gene.annot.xls  ';
			}
			if ($option{$cluster}{blast_option}=~/-kegg/) {
				print EXP ' -kegg '.$blast_xls_fa_name.'.ko  ';
			}
			if ($option{$cluster}{blast_option}=~/-nt/){
				print EXP ' -nt '.$blast_xls_fa_name.'.blast.Nt.xls  ';
			}
			print EXP ' && echo finish add blast annot to '.$exp_flag.' at time '.$time."\n";
			if (exists $option{$cluster}{nr2go_option}) {
				print EXP 'perl '.$config{"annot2goa"}.'   '.$godir.'/'.$gene.'.blast.Nr.xml.annot  '.$diff_dir.'/'.$gene_name.'  && echo finish annot2go at time '.$time."\n";
				print EXP 'perl '.$config{"addGOAnnot"}.'   -go '.$diff_dir.'/'.$gene_name.' -input  '.$diff_dir.'/'.$exp_flag.'_annot.xls -output '.$diff_dir.'/annotation.xls   && echo finish add go to '.$exp_flag.' at time '.$time."\n";
				print EXP 'perl '.$config{"annot_version"}.' '.$lib.' '.$diff_dir.'/'.'database.version.txt'."\n";
			}
		}
		close EXP;

		my $annotation_shell = &_Upload("$diff_dir/annotation.xls","$up_blast_dir/annotation.xls");
		print UPLOADSH $annotation_shell;
	if ((exists $option{$cluster}{Diff}) && (scalar(@samplename)>1) ) {
		my $up_diff_dir=$upload_dir.'/diff_exp';
		system "mkdir -p  $up_diff_dir" unless(-d $up_diff_dir);
		system "mkdir -p  $up_diff_dir/GO" unless (-d "$up_diff_dir/GO");
		system "mkdir -p  $up_diff_dir/GO/Category "unless (-d "$up_diff_dir/GO/Category");
		system "mkdir -p  $diff_dir/pathway " unless (-d "$diff_dir/pathway");
		system "mkdir -p  $up_diff_dir/Pathway " unless (-d "$up_diff_dir/Pathway");
		system "mkdir -p  $up_diff_dir/Genelist" unless (-d "$up_diff_dir/Genelist"); 
		my $pvaluedir=$shdir.'/PvalueSH';##by dushuo 2010-12-9
		system "mkdir -p  $pvaluedir" unless (-d $pvaluedir);##by dushuo 2010-12-9

	  ############# pvalue
		my $glist='';
		my @vs=split /\,/, $option{$cluster}{Diff};
		open TOTALSH,'>',$pvalue_sh or die "can't open the pvalue sh $pvalue_sh";##by dushuo 2010-12-8

		my $komap=$config{map};
		if ($option{$cluster}{blast_option}=~/-kegg\s*(\S+)/) {
			my $tmp='map_'.$1;
			$komap=$config{$tmp};
		}

		foreach my $vspart (@vs) {
			my ($nameA,$nameB)=split /&/,$vspart;
			my $vsname = $nameA.'-vs-'.$nameB;
			my ($numA, $numB);
			for(my $tmp_num = 0; $tmp_num < @samplename; $tmp_num++){
				if($nameA eq $samplename[$tmp_num]){
					$numA = $tmp_num + 1;
				}
				if($nameB eq $samplename[$tmp_num]){
					$numB = $tmp_num + 1;
				}
			}
			my $vsdir=$pvaluedir.'/'.$vsname;
			system "mkdir -p  $vsdir" unless (-d $vsdir);##by dushuo 2010-12-9
			print TOTALSH "sh $pvaluedir/$vsname/$vsname.sh\n";##by du
			open PVALUESH,'>',"$pvaluedir/$vsname/$vsname.sh";##by du
			print PVALUESH "export LD_LIBRARY_PATH=$rna_lib:\$LD_LIBRARY_PATH\n";#add by tangqq
			my $tmp_alg = $exp_flag.'2Pvalue';
			print PVALUESH 'perl  '.$config{$tmp_alg};
			print PVALUESH  '  '.$diff_dir.'/'.$exp_flag.'.xls   '.$numA.'   '.$numB.'  '.$diff_dir.'/'.$vsname.'.xls    && echo finish get '.$vsname.'  pvalue  at time '.$time."\n";
			print PVALUESH 'perl  '.$config{FDR}.'   '.$diff_dir.'/'.$vsname.'.xls  '.$diff_dir.'/'.$vsname.'.compare.xls  && echo finish '.$vsname.' FDR at time '.$time."\n";
			print PVALUESH  'perl '.$config{get_diff_from_pvalue}.' -pvalue '.$diff_dir.'/'.$vsname.'.compare.xls  -go '.$godir.'/'.$gene.'.blast.Nr.xml.wego  -outdir '.$diff_dir.'     -ko  '.$blast_xls_fa_name.'.ko -path1 '.$path1.' -path2 '.$path2.' -shdir '.$vsdir.' && echo finish get '.$vsname.' diff go and ko at time '.$time."\n";
			print PVALUESH 'awk '."'{print ".'$1"\t"$7}'."' ".$diff_dir.'/'.$vsname.'.diff.xls | sed '."'".'1c Gene\t'.$vsname."'".' > '.$diff_dir.'/'.$vsname.'.glist   && echo finish get '.$vsname.'.glist at time  '.$time."\n";
			print PVALUESH 'perl '.$config{plotDiffExp}.' -input '.$diff_dir.'/'.$vsname.'.compare.xls  -Head  -output  '.$diff_dir.'/'.$vsname.'.svg  -nameA  '.$samplename[$numA-1].' -nameB '.$samplename[$numB-1].' -fdr  0.001  -log2  1 -exp_alg '.$exp_flag.' && echo finish draw diffExp '.$vsname.'  at  time '.$time."\n";
			print PVALUESH <<SHELL;
if [ -L "$diff_dir/pathway/" ]; then rm -rf "$diff_dir/pathway/*"; else mkdir -p $diff_dir/pathway/; fi
SHELL
			print PVALUESH 'perl '.$config{pathfind}.' -fg '.$diff_dir.'/'.$vsname.'.ko -bg '.$blast_xls_fa_name.'.ko -output '.$diff_dir.'/pathway/'.$vsname.'.path  -komap '.$komap.'   && echo finish '.$vsname.' pathfind at time '.$time."\n";
			print PVALUESH 'awk '."'{print ".'$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t"$10 }'."' ".$diff_dir.'/'.$vsname.'.diff.xls  >  '.$diff_dir.'/pathway/'.$vsname.'.diff.xls  && echo finish get pathway '.$vsname.'.diff.xls at time '.$time."\n";
			print PVALUESH <<MAPSHELL;
if [ -d "$diff_dir/pathway/$vsname\_map" ]; then rm -rf $diff_dir/pathway/$vsname\_map; fi
mkdir -p $diff_dir/pathway/$vsname\_map
MAPSHELL
			print PVALUESH 'perl '.$config{keggMap}.' -ko '.$diff_dir.'/'.$vsname.'.ko -diff '.$diff_dir.'/pathway/'.$vsname.'.diff.xls -outdir '.$diff_dir.'/pathway/'.$vsname.'_map  -komap '.$komap.' && echo finish  '.$vsname.' keggmap at time '.$time."\n";
			print PVALUESH 'cp '.$diff_dir.'/'.$vsname.'.ko '.$diff_dir.'/pathway/'."\n";
			my @diffs;
			push @diffs,"$diff_dir/$vsname.compare.xls","$diff_dir/$vsname.diff_annot.xls","$diff_dir/diff_gene_counts.png","$diff_dir/$vsname.svg","$diff_dir/$vsname.png";
			foreach my $diff_result(@diffs){
				my $diff_filename = basename($diff_result);
				my $diff_target = $up_diff_dir.'/Genelist/'.$diff_filename;
				my $diff_script = &_Upload($diff_result,$diff_target);
				print UPLOADSH $diff_script;
			}
			my @go_category;
			push @go_category,"$diff_dir/annotation/$vsname.GO.svg","$diff_dir/annotation/$vsname.GO.png";
			foreach my $go_ca_result(@go_category){
				my $go_ca_filename = basename($go_ca_result);
				my $go_ca_target = $up_diff_dir.'/GO/Category/'.$go_ca_filename;
				my $go_ca_script = &_Upload($go_ca_result,$go_ca_target);
				print UPLOADSH $go_ca_script;
			}
			$glist.=$diff_dir.'/'.$vsname.'.glist,';
			print UPLOADSH 'cat '.$godir.'/wego.head  '.$diff_dir.'/annotation/'.$vsname.'.'.$gene.'.blast.Nr.xml.wego  > '.$up_diff_dir.'/GO/Category/'.$vsname.'.gene2GO.xls  && echo finish upload '.$vsname.'.gene2Go.xls at time '.$time."\n";
			my $go_cat_shell = &_Upload("$diff_dir/annotation/$vsname.GO.xls","$up_diff_dir/GO/Category/$vsname.GO2gene.xls");
			print UPLOADSH $go_cat_shell;
			close PVALUESH;##by du
		}
		close TOTALSH;
		system "mkdir -p  $diff_dir/diff_go" unless(-d "$diff_dir/diff_go");
		$glist=~s/,$//;
		open FUNCTIONALSH,'>',$functional_sh or die "can't open the functional sh  $functional_sh";
		print FUNCTIONALSH <<FUNC;
if [ -d "$diff_dir/diff_go/" ]; then rm -rf $diff_dir/diff_go/GO/*; fi
FUNC
		print FUNCTIONALSH 'perl '.$config{functional}.' -go  -glist '.$glist.'  -sdir '.$diff_dir.' -species '.$gene_name.' -outdir '.$diff_dir.'/diff_go  && echo finish  functional diffgo at time '.$time."\n";
		print FUNCTIONALSH 'rm -rf  '.$diff_dir.'/diff_go/GO/plugins/.svn  && echo finish rm plugins svn at time '.$time."\n";
		print FUNCTIONALSH 'perl '.$config{genPathHTML}.' -indir '.$diff_dir.'/pathway   && echo finish diffpathway at time '.$time."\n";
		print FUNCTIONALSH 'perl '.$config{diff_exp_annot}.' '.$diff_dir.'  && echo finish add annotation infomation to DEGs at time '.$time."\n";
		close FUNCTIONALSH;
		my $diff_gos_shell = &_Upload("$diff_dir/diff_go/GO","$up_diff_dir/GO/EnrichmentAnalysis");
		print UPLOADSH $diff_gos_shell;
		print CLEARSH 'rm -rf '.$diff_dir.'/diff_go/GO   && echo finish rm diffgo\(dge\s\) mid file at time '.$time."\n";

		close CLEARSH;
		my $diff_path_shell = &_Upload("$diff_dir/pathway","$up_diff_dir/Pathway");
		print UPLOADSH $diff_path_shell;
	}
}


##########################################################################################
################################# SSR Analysis ###########################################
##########################################################################################
my $ssr_sh = $shdir.'/ssr.sh';
if(exists $option{$cluster}{'SSR'}){
	my $ssr_dir = "$outdir/SSR";
	system("mkdir -p $ssr_dir");
	open SSR,">$ssr_sh" or die "Can't open file $ssr_sh!";
	print SSR<<SSRTEXT;
rm -rf $ssr_dir/*
cd $ssr_dir
ln -s $clu_dir/$gene $ssr_dir/$gene
perl $config{misa} $ssr_dir/$gene $ssr_dir $gene_name $option{$cluster}{SSR}
perl $config{'2primer_designer'} $ssr_dir/$gene_name.ssr.out.xls $ssr_dir/$gene_name.raw_primer.out.xls $ssr_dir/$gene_name.primer_results.out.xls
perl $config{'3filter_primer'} $ssr_dir/$gene_name.primer_results.out.xls $ssr_dir/$gene_name.rescreen.out.xls $ssr_dir/$gene_name.blastin
/opt/blc/genome/biosoft/blast-2.2.23/bin/formatdb -i $ssr_dir/$gene -p F
/opt/blc/genome/biosoft/blast-2.2.23/bin/blastall -i $ssr_dir/$gene_name.blastin  -d $ssr_dir/$gene -p blastn -o $ssr_dir/$gene_name.blast.out.xls -F F -a 15
perl $config{'5blast_parse'} $ssr_dir/$gene_name.blast.out.xls $ssr_dir/$gene_name.query_sbjct.out.xls $ssr_dir/$gene_name.stactic.out.xls 4
perl $config{'6generate_primer_pair_new'} $ssr_dir/$gene_name.query_sbjct.out.xls $ssr_dir/$gene_name.primer.tab.xls
perl  $config{'7final_pair'} $ssr_dir/$gene_name.primer.tab.xls  $ssr_dir/$gene_name.rescreen.out.xls  $ssr_dir/$gene_name.inter_primer.tab.xls
perl $config{'8product_ssr_check'}  $ssr_dir/$gene  $ssr_dir/$gene_name.inter_primer.tab.xls  $ssr_dir/$gene_name.product_file.out.xls
perl $config{ssr_finder} --ssr 12 $ssr_dir/$gene_name.product_file.out.xls  $ssr_dir/$gene_name.product_ssr.out.xls
perl $config{'9filter_ssr'}  $ssr_dir/$gene_name.product_ssr.out.xls  $ssr_dir/$gene_name.inter_primer.tab.xls  $ssr_dir/$gene_name.only_primer.tab.xls
perl $config{'10get_fitprimer'} $ssr_dir/$gene_name.only_primer.tab.xls  $ssr_dir/$gene_name.product_ssr.out.xls  $ssr_dir/$gene_name.rescreen.out.xls  $ssr_dir/$gene_name.final_primer.out
sort -k 1,1 $ssr_dir/$gene_name.final_primer.out > $ssr_dir/$gene_name.final_primer.out.xls

echo Final primered SSR:
awk '{print \$1}' $ssr_dir/$gene_name.final_primer.out.xls | awk 'BEGIN{FS="_"}{print \$1}' | sort -u | wc -l
echo Final product:
wc -l $ssr_dir/$gene_name.final_primer.out.xls

perl $config{stat_draw} $ssr_dir/$gene_name.statistics.xls $ssr_dir
SSRTEXT
	close SSR;
	my $up_ssr_dir = "$upload_dir/SSR";
	system("mkdir -p $up_ssr_dir") unless(-d $up_ssr_dir);
	my @ssr_arr;
	push @ssr_arr, "$ssr_dir/$gene_name.misa.xls", "$ssr_dir/$gene_name.statistics.xls", "$ssr_dir/$gene_name.primer_results.out.xls", "$ssr_dir/$gene_name.final_primer.out.xls", "$ssr_dir/SSR_statistics.png", "$ssr_dir/primerin_example.txt", "$ssr_dir/ssr_stat.xls";
	for my $ssr_result(@ssr_arr){
		my $ssr_filename = basename($ssr_result);
		my $ssr_target = $up_ssr_dir."/".$ssr_filename;
		my $ssr_shell = &_Upload($ssr_result, $ssr_target);
		print UPLOADSH $ssr_shell;
	}
	undef @ssr_arr;
}

##########################################################################################
################################# SNP Analysis ###########################################
##########################################################################################
#my $snp_sh = $shdir.'/snp.sh';
my $snp_pl = $shdir.'/soapsnp.pl';
if(exists $option{$cluster}{'SNP'}){
	my $snp_dir = "$outdir/SNP";
	system("mkdir -p $snp_dir");
	my $up_snp_dir = "$upload_dir/SNP";
	system("mkdir -p $up_snp_dir");
	my $snp_sh_dir = $shdir.'/SNPSH';
	system("mkdir -p $snp_sh_dir");

	my $sam_snp_join = "";
	open SNP_PL,">$snp_pl";
	print SNP_PL<<SNPTEXT;
#!/usr/bin/perl
use strict;
use lib "/opt/rocks/lib/perl5/5.10.1";
use Thread 'async';

SNPTEXT

	my $snp_file_list = "";
	my $cns_dir_list = "";
	my $snp_sam_num = scalar @samplename;
	for my $sam_snp_name(@samplename){
		my $sam_snp_sh = $snp_sh_dir.'/'.$sam_snp_name.'.snp.sh';
		my $sam_snp_dir = "$snp_dir/$sam_snp_name";
		system("mkdir -p $sam_snp_dir");
		my $max_read_length = $option{$sam_snp_name}{length};
		my $finalfilter = '-z gz -m 15 -n 2 -q 20 -l '.$max_read_length;
#		$snp_file_list .= ",$up_snp_dir/$sam_snp.snp.xls";
#		$snp_file_list .= ",$sam_snp_dir/snp_filter/$sam_snp_name.snp.xls";
		$snp_file_list .= ",$sam_snp_dir/snp_filter/$sam_snp_name.All_chr.snp.xls";
		$cns_dir_list .= ",$sam_snp_dir/cns_filter";
#		print "$config{individual_snp_run}\n";
		system("perl $config{individual_snp_run} -op_snp \'$option{$cluster}{'SNP'} -L $max_read_length\' -finalfilter \'$finalfilter\' -queue $queue1 -Pname $option{$cluster}{sub_prjct_id} -fa $clu_dir/$gene -soap $outdir/soap/$sam_snp_name.PESoap.gz,$outdir/soap/$sam_snp_name.PESoapSingle.gz -name $sam_snp_name -outdir $snp_dir/$sam_snp_name -final $up_snp_dir -num_fa 50000 -sam_num $snp_sam_num &>/dev/null");
		open SAMSNP,">$sam_snp_sh";
		print SAMSNP "sh $sam_snp_dir/SH/individual_snp.$sam_snp_name.sh\n";
#		if(scalar(@samplename) == 1){
#			print SAMSNP "rm -rf $sam_snp_dir/cns_filter";
#		}
		close SAMSNP;
		my $sam_snp = $sam_snp_name;
		$sam_snp =~ s/\W+/\_/g;
		$sam_snp =~ s/^(\d+)/S_$1/g;
		print SNP_PL<<SUBSNP;
my \$$sam_snp = async {
	my \$step1=1;
	\$step1=system ("sh $sam_snp_sh");
	if(\$step1 eq 0){return 0}else{return 1}
};

SUBSNP
		$sam_snp_join .= "\$$sam_snp->join() == 0 && ";
	}
	$sam_snp_join =~ s/ \&\& $//g;
	print SNP_PL<<SNPTEXT2;
if ($sam_snp_join){
	print "finish sample soapsnp!";
	exit 0;
}
else{
	print "Warn soapsnp wrong!";
	exit 1;
}
SNPTEXT2
	close SNP_PL;
	$snp_file_list =~ s/^,//g;
	$cns_dir_list =~ s/^,//g;
	if(@samplename > 1){
		open SNPDIFF,">$shdir/soapsnp_in_group.sh";
		print SNPDIFF "sh $config{snp_diff} $up_snp_dir/All-Sample $snp_file_list $cns_dir_list ".scalar(@samplename)."\n";
		$cns_dir_list =~ s/,/ /g;
		print SNPDIFF "rm -rf $cns_dir_list\n";
		close SNPDIFF;
	}
}

###########################################################################################
###############  生成 结题 报告
###########################################################################################
system "mkdir -p  $upload_dir/image" unless (-d "$upload_dir/image");

print UPLOADSH 'cp  '.$config{"image"}.'/*  '.$upload_dir.'/image/  && echo finish upload html image at time '.$time."\n";
#print UPLOADSH 'perl  '.$config{denovo_readme}.' -lib  '.$lib.' -dir '.$upload_dir.'    && echo finish get readme at time '.$time."\n";
#print UPLOADSH 'perl  '.$config{denovo_readme_en}.' -lib  '.$lib.' -dir '.$upload_dir.'    && echo finish get readme at time '.$time."\n";
#print UPLOADSH 'perl  '.$config{check}." -fdir $upload_dir/image -ftype gif,jpg,png -n 1,13,1  -out upload.e.info \n";
print UPLOADSH 'perl  /ifs4/BC_PUB/biosoft/pipe/bc_ba/RNAdenovo/Denovo_Pipeline.2.0/genHTML4Denovo_2.0_for_zhou.pl -indir '.$upload_dir.' -lib '.$lib.' -outdir '.$upload_dir.'    && echo finish get report html at time '.$time."\n";
print UPLOADSH 'perl  '.$config{svn_version}.' /ifs4/BC_PUB/biosoft/pipe/bc_ba/RNAdenovo/Denovo_Pipeline.2.0/genHTML4Denovo_2.0_for_zhou.pl '.$config{image}.' '.$config{image}.'  '.$config{denovo_quality}.' '.$config{fq_num_bp}.' '.$config{get_chosen_fa}.' '.$config{denovo_readme}.' '.$config{denovo_readme_en}.' > '.$version_dir.'/upload.txt'."\n";# chuanjun
print UPLOADSH "sh $outdir/upload_tar.sh\n";
print UPLOADSH "sh $outdir/SH/rna_data.sh\n";
close UPLOADSH;
##########################################################################################
############# qsub
##########################################################################################
my $step=1;
my $final_sh=$outdir.'/'.$cluster.'.final.sh';
open FINALSH,'>',$final_sh or die "can't open the sh $final_sh";
print FINALSH 'echo stat at '.$time."\n\n";
$cwdmv='find  '.$shdir.' -name "*.sh.*" -maxdepth 1  -exec mv {} '.$oldshdir.' ";"';
print FINALSH $cwdmv."\n\n";

#########added by dushuo 11-22-2010
print FINALSH '######## '.$step.',sample assembly ( maybe max when as high as 20G )'."\n";
open ASYNC,'>',"$shdir/sample_assembly.sh";
for my $sam(@samplename){
	my $sam_sh="$samsh/$sam.sample_assembly.sh";
	print ASYNC "sh $sam_sh\n";
}
close ASYNC;
$step++;
my $asm_memory;
if(exists $option{$cluster}{'CPU'}){
	$asm_memory = $option{$cluster}{'CPU'} * 3 + 5;
}
else{
	$asm_memory = 7 * 3 + 5;
}
print FINALSH 'perl  '.$config{qsub_sge}.' --getmem  --jobprefix sample_assembly --queue '.$queue1.' --subprjctid '.$option{$cluster}{sub_prjct_id}.' --resource vf=3G --lines  1  --convert no ',$shdir,"/sample_assembly.sh && echo finish filter dup as f1 sam_clu at $time\n\n";
######### check assembly result ###########################
print FINALSH<<CHECKAS;
if [ "\$(cat \$(ls ./SH/sample_assembly.sh.*.log| awk '\$0!~/mem.log/'))" = "All jobs finished!" ];then
	echo There is no error in assembly!
else
	echo Some error was happened in assembly!
	exit 1
fi
CHECKAS
print FINALSH 'cd ',$shdir,' && qsub -cwd -S /bin/bash -q '.$queue1.' -P '.$option{$cluster}{sub_prjct_id}.' -l vf=300M '.$quality_sam_sh.' -e '.$quality_sam_sh.'.e -o '.$quality_sam_sh.'.o'."\n\n";
print VERSION 'perl  '.$config{svn_version}.'  '.$config{denovo_quality}.' > '.$version_dir.'/sam_quality.txt '."\n";
###################
if (-s $clustersh ) {
	print FINALSH '######## '.$step.',  cluster  ( use tgicl  ; maybe need 1G ; ) '."\n";
	##by chuanjun
	print VERSION 'perl  '.$config{svn_version}.'      '.$config{get_chosen_fa}.'  '.$config{replace_geneid}.'  '.$config{fa_quality}.' '.$config{barplot}.' '.$config{qsub_sge}.'  > '.$version_dir.'/cluster.txt '."\n";
	print FINALSH 'perl  '.$config{qsub_sge}.'  --reqsub  --getmem  --jobprefix clu --queue '.$queue1.' --subprjctid '.$option{$cluster}{sub_prjct_id}.' --resource vf=4G --lines  100  --convert no '.$clustersh.' &&  echo finish cluster  clu sh  at  '.$time."\n\n";
	print FINALSH 'cd ',$shdir,' && qsub -cwd -S /bin/bash -q '.$queue1.' -P '.$option{$cluster}{sub_prjct_id}.' -l vf=300M '.$quslity_sh.' -e '.$quslity_sh.'.e -o '.$quslity_sh.'.o'."\n\n";
	$step++;
}
############## for paralleling ######all things about ANA is for paralleling -->dushuo

print FINALSH '######## '.$step.', annot  '."\n";
print FINALSH "perl $shdir/annot.pl\n\n";
$step++;

my $exp_cov="$shdir/$exp_flag\_cov.pl";
my $anacmd='';
open ANA,'>',$exp_cov;
print ANA<<ALYSIS;
#!/usr/bin/perl
use strict;
use lib "/opt/rocks/lib/perl5/5.10.1";
use Thread 'async';

ALYSIS

print FINALSH '######## '.$step.', '.$exp_flag.' && coverage  '."\n";
print FINALSH "perl $shdir/$exp_flag\_cov.pl\n\n";
$step++;

my $annot_analysis = "";
if (exists $option{$cluster}{blast_option}) {
	$annot_analysis .= " && \$blast_go->join() == 0";
###############--> dushuo 2010-11-24
	print ANNO<<LYSIS;
	my \$step2=1;
	\$step2=system ("perl $shdir/cds_go.pl") if (\$step1 eq 0);
LYSIS
	open CDSGO,'>',"$shdir/cds_go.pl";
	print CDSGO<<GOCDS;
#!/usr/bin/perl
use strict;
use lib "/opt/rocks/lib/perl5/5.10.1";
use Thread 'async';

GOCDS
	my $judge='';
	print VERSION 'perl  '.$config{svn_version}.'  '.$config{search_database}.'  '.$config{finish}.'    '.$config{cog_R}.'  '.$config{fastaDeal}.'   '.$config{blast_parser}.' '.$config{cog_parser}.'   '.$config{blast2ko}.' '.$config{pathfind}.' '.$config{keggMap_nodiff}.' '.$config{genPathHTML}.' > '.$version_dir.'/database_annot.txt '."\n";
	if (exists $option{$cluster}{nr2go_option}) {
		print VERSION 'perl  '.$config{svn_version}.'  '.$config{blast_m0_m7_zhd}.'  '.$config{annot2wego}.'  '.$config{drawGO}.'    '.$config{qsub_sge}.'    > '.$version_dir.'/go.txt '."\n";
		$judge.="\$go->join() == 0 &&";
		print CDSGO 'my $go = async {'."\n\tmy \$step1=1;\n".'	$step1=system ("perl '.$config{qsub_sge}.' --reqsub --getmem --queue '.$queue1.' --subprjctid '.$option{$cluster}{sub_prjct_id}.' --resource vf=15g --jobprefix blast2go --lines 10 --interval 300 '.$shdir.'/'.$gene.'.blast2go.sh  && echo finish blast2go at '.$time.'");'."\n";
		print CDSGO '	my $step2=1;'."\n".'	$step2=system ("perl  '.$config{qsub_sge}.' --convert no --reqsub --getmem  --jobprefix go --queue '.$queue1.' --subprjctid '.$option{$cluster}{sub_prjct_id}.' --resource vf=6G --lines  100   '.$go_sh.'  &&  echo finish go at  '.$time."\") if(\$step1 eq 0);\n";
		print CDSGO "\tif(\$step2 eq 0){return 0}else{return 1}\n};\n";
	}
	print VERSION 'perl  '.$config{svn_version}.'    '.$config{"get_cds_blast.pl"}.'    '.$config{prepare_data}.'  '.$config{build_model}.'  '.$config{estscan}.'  '.$config{fa_quality}.' '.$config{barplot}.' > '.$version_dir.'/CDS.txt '."\n";

	print CDSGO 'my $cds = async {'."\n\tmy \$step1=1;\n".'	$step1=system ("perl  '.$config{qsub_sge}.' --reqsub --getmem --jobprefix cds --queue '.$queue1.' --subprjctid '.$option{$cluster}{sub_prjct_id}.' --resource vf=0.5G --lines  100   '.$cds_sh.' &&  echo finish CDS at  '.$time."\");\n";
	print VERSION 'perl  '.$config{svn_version}.'  '.$config{sign}.'  '.'  > '.$version_dir.'/sign.txt '."\n";
	print CDSGO "\tmy \$step2=1;\n".'	$step2=system ("sh '.$fs_sign_sh.' &&  echo finish gene orientation at  '.$time."\") if (\$step1 eq 0);\n\tif(\$step2 eq 0){return 0}else{return 1}\n};\n\n";
	$judge.=" \$cds->join() == 0";
	print CDSGO<<GOCDS;
if($judge){
	print "Finish cds go  !";
	exit 0;
}else{
	print "Warn go || cds || sign wrong!";
	exit 1;
}
GOCDS
	close CDSGO;
}

print ANNO<<LYSIS;
	if(\$step2 eq 0){return 0}else{return 1}
};

LYSIS
if (exists $option{$cluster}{soap}) {
	$annot_analysis .= " && \$soap->join() == 0";
	print VERSION 'perl  '.$config{svn_version}.'  '.$config{fa_quality}.'  '.$config{ReadsRandomInGene}.'  '.$config{Soap_Coverage}.' '.$config{qsub_sge}.' > '.$version_dir.'/soap.txt '."\n";
	
	print ANNO 'my $soap = async {'."\n\tmy \$step1=1;\n".'	$step1=system ("perl  '.$config{qsub_sge}.' --reqsub --getmem --jobprefix bwt --queue '.$queue1.' --subprjctid '.$option{$cluster}{sub_prjct_id}.' --resource vf=2G --lines  100   '.$bwt_sh.' &&  echo finish 2bwt at  '.$time."\");\n";
	print ANNO "\tmy \$step2=1;\n".'	$step2=system ("perl  '.$config{qsub_sge}.' --reqsub --getmem  --jobprefix soap --queue '.$queue1.' --subprjctid '.$option{$cluster}{sub_prjct_id}.' --resource vf=2G --lines  1   '.$soap_sh.' &&  echo finish soap at  '.$time."\") if (\$step1 eq 0);\n";
	if(exists $option{$cluster}{SNP}) {
		my $snp_result_files = "";
		foreach (@samplename){
			$snp_result_files .= ",$upload_dir/SNP/$_.snp.xls";
		}
		$snp_result_files =~ s/^,//g;
		if(@samplename == 1){
			print ANNO "\tmy \$step3=1;\n\t\$step3=system (\"perl $shdir/soapsnp.pl\") if (\$step2 eq 0);\n";
			print ANNO "\tmy \$step4=1;\n\t\$step4=system (\"perl $config{snp_statistics} -snp $snp_result_files -out $upload_dir/SNP && echo finish statistics snp at $time\") if (\$step3 eq 0);\n";
			print ANNO "\tif(\$step4 eq 0){return 0}else{return 1}\n};\n\n";
		}
		elsif(@samplename > 1){
			print ANNO "\tmy \$step3=1;\n\t\$step3=system (\"perl $shdir/soapsnp.pl\") if (\$step2 eq 0);\n";
			print ANNO "\tmy \$step4=1;\n\t\$step4=system (\"perl $config{qsub_sge} --getmem --jobprefix soapsnp.diff --queue $queue1 --subprjctid $option{$cluster}{sub_prjct_id} --resource vf=2G --lines 100 $shdir/soapsnp_in_group.sh &&  echo finish soapsnp diff analysis at  $time\") if (\$step3 eq 0);\n";
			print ANNO "\tmy \$step5=1;\n\t\$step5=system (\"perl $config{snp_statistics} -snp $snp_result_files -out $upload_dir/SNP && echo finish statistics snp at $time\") if (\$step4 eq 0);\n";
			print ANNO "\tif(\$step5 eq 0){return 0}else{return 1}\n};\n\n";
		}
	}
	else{
		print ANNO "\tif(\$step2 eq 0){return 0}else{return 1}\n};\n\n";
	}
	print ANA "my \$coverage = async {\n\tmy \$step1=1;";
	print ANA "\t\$step1=system(\"perl $config{qsub_sge} --reqsub --getmem --jobprefix cover --queue $queue1 --subprjctid $option{$cluster}{sub_prjct_id} --resource vf=2G --lines  1  $coverage_sh &&  echo finish soap coverage at  \\`date +\%y-\%m-\%d.\%H:\%M:\%S\\`\");\n\tif(\$step1 eq 0){return 0}else{return 1}\n};\n\n";

	print ANA "my \$analysis = async {\n";
	print VERSION 'perl  '.$config{svn_version}.'  '.$config{$exp_flag}.'  '.$config{get_sam_by_exp_new}.' '.$config{annot2goa}.'  '.$config{addGOAnnot}.'    > '.$version_dir.'/'.$exp_flag.'.txt '."\n";

	$anacmd.="\tmy \$step1=1;\n\t\$step1=system(\"perl  $config{qsub_sge} --reqsub --getmem  --jobprefix $exp_flag --queue $queue2 --subprjctid $option{$cluster}{sub_prjct_id}  --resource vf=15G --lines  100 $exp_sh  &&  echo finish $exp_flag at  \\`date +\%y-\%m-\%d.\%H:\%M:\%S\\`\");\n";
	if (defined $option{$cluster}{Diff}) {
		print VERSION 'perl  '.$config{svn_version}.'  '.$config{exp2Pvalue}.'  '.$config{FDR}.' '.$config{get_diff_from_pvalue}.'  '.$config{plotDiffExp}.'  '.$config{pathfind}.'  '.$config{keggMap}.'  '.$config{functional}.'  '.$config{genPathHTML}.'  > '.$version_dir.'/pvalue.txt '."\n";
		$anacmd.="\tmy \$step2=1;\n\t\$step2=system(\"rm -rf $outdir/diff/pathway/*\") if (\$step1 eq 0);\n";
		$anacmd.="\tmy \$step3=1;\n\t\$step3=system(\"sh $pvalue_sh  &&  echo finish pvlaue at \\`date +\%y-\%m-\%d.\%H:\%M:\%S\\`\") if(\$step2 eq 0);\n";##changed by dushuo from line 10 to 1
		$anacmd.="\tmy \$step4=1;\n\t\$step4=system(\"perl  $config{qsub_sge} --reqsub --getmem   --convert no --jobprefix functional --queue $queue2 --subprjctid $option{$cluster}{sub_prjct_id}  --resource vf=2G --lines  200 $functional_sh  &&  echo finish functional  at  \\`date +\%y-\%m-\%d.\%H:\%M:\%S\\`\") if (\$step3 eq 0);\n";
		$anacmd.="\tif(\$step4 eq 0){return 0}else{return 1}\n};\n\n";
	}else{
		$anacmd.="\tif(\$step1 eq 0){return 0}esle{return 1}\n};\n\n";
	}
	print ANA $anacmd;
	print ANA<<LYSIS;
if ( \$coverage->join() == 0 && \$analysis->join() == 0 ) {
	print "analysis of All-Unigene done!\\n";
	exit 0;
}else{
	print "Warn coverage or analysis may end in error";
	exit 1;
}
LYSIS
}

if(exists $option{$cluster}{'SSR'}){
	$annot_analysis .= " && \$ssr->join() == 0";
	print ANNO<<ONNA;
my \$ssr = async {
	my \$step1=1;
	\$step1=system (\"perl $config{qsub_sge} --reqsub --getmem --jobprefix ssr --queue $queue1 --subprjctid $option{$cluster}{sub_prjct_id} --resource vf=4G --lines  100 $ssr_sh &&  echo finish ssr at  \`date +\%y-\%m-\%d.\%H:\%M:\%S\`\");
	if(\$step1 eq 0){return 0}else{return 1}
};

ONNA
}

$annot_analysis =~ s/^ && //g;
print ANNO<<ONNA;
if ($annot_analysis){
	print "finish database, soap and ssr!";
	exit 0;
}else{
	print "Warn database, soap or ssr wrong!";
	exit 1;
}
ONNA
print FINALSH '######## '.$step.',  upload   (maybe need 300M   ; ) '."\n";
print FINALSH 'cd ',$shdir,'  && qsub -cwd -S /bin/bash -q ',$queue1 ,' -P ',$option{$cluster}{sub_prjct_id} ,' -l vf=500M ',$upload_sh,"\n";
print FINALSH 'echo  all finish  at  '.$time."\n";
close FINALSH;

open TARSH,'>',$outdir.'/upload_tar.sh' or die "can't open the sh of tar upload $outdir/upload_tar.sh ";
print TARSH "cd $outdir \n";
print TARSH "/bin/tar -cjvf $cluster.tar.bz2  $cluster \n";
print TARSH $config{md5sum}."  $cluster.tar.bz2 > $cluster.tar.bz2.md5 \n";
close TARSH;
my $usage=<<USAGE;
please  follow this step:

1, start && you can choose from these two ways to run your process
cd $outdir

nohup sh $cluster.final.sh 1> $cluster.final.sh.nohup 2> $cluster.final.sh.error &

2, Clean up mid files to free diskspace
if you are sure the all finish and want to free diskspace ,then

cd $shdir
qsub -cwd -S /bin/bash -q $queue1 -l vf=300M clear.sh 

all finish
USAGE
print STDOUT $usage;
open HELP,'>',$outdir.'/help.txt' or die "can't open the $outdir/help.txt ";
print HELP $usage;
close HELP;

###
sub fq_suffix{
	my ($sam,$option)=@_;
	my $suffix="";
	if ($option=~/\s(-t|--type)\s+(\d)/) {
		if ($2==1) {
			$suffix.=".gz";
		}elsif ($2==2){
			$suffix.=".bz2";
			print STDERR "error: soap2.21 unrecognized bz2 file \n \ttest dir: /ifs1/DGE_SR/nixiaoming/work/Find_bug_denovo/2010-10-26.filter_fq/test/soap.bz2\n";
		}elsif ($2==0){
		}else{
			print STDERR "error: $sam filter_fq_option -t \n";
		}
	}else{
		$suffix.=".gz";
	}
	return $suffix;
}
#######################################
# to return the words to judge the file is exists or not
# we move the result file to the upload directory, and then make a link(ln -s) to the old place, so we should check the file is exists or not
# if exists, do rm, or skip
#######################################
sub _Upload{
	my $source = shift;  # the result
	my $target = shift;  # the upload file
	my $shell = <<TEXT;
if [ -L "$source" ]; then
	echo The file or directory: $source is a link
elif [ -e "$source" ]; then
	if [ -e "$target" ]; then rm -rf $target; fi && mv $source $target && ln -s $target $source
else
	echo ERROR:  The file $source is not exists! The move operation failed! Please check it!
fi
TEXT
	return $shell;
}


my $numlib="$outdir\/numlib";
open NUMLIB,'>',$numlib;
print NUMLIB <<NUM;
######################################################################
################ how to configure this lib to fit your need ##########
######################################################################

#1.here is a example of lib , you can modify it as your need
#2.default not to stop or wait , and if you want ,you can defined 'stop=stop' or 'stop=wait' in a new line in the step area 
#3.I don't encourage you to change the 'depth or xdepth or fdir' when the lib now is well used 

######################## what you can change #########################
#1.if you don't want to check a certain step , then annotate all lines about the step or delete them 
#2.if you don't want to check one certain type of files , then remove it and remove its corresponding 'n'(file number) as well 

########you may be confused about the 'n' with 's|v|nothing' #########
#1.if the number of file belongs to sample number , then you can write the files num as 1s , s represents 'Sample'
#2.if the number of file belongs to the vs number , then you can write the files num as 1v , v represents 'vs'
#3.the file number can be result of calculation , I can do maths like '2s+2'|'2+2s'|'2s+3v'|'3v+2'|'2+3v'|'3v+2s'|'2s'|'3v'|'5'

# never add annotation at the end of one line , I can't recognise it #

###################these steps contains more steps , we can write them as this form
>sample_assembly
sample_assembly=filter,dup,as,f1,sam_clu

>cluster
fdir=cluster
ftype=R,Z,cidx,cluster,dat,fa,fa_cl_clusters,gnu,list,log,lst,nhr,nin,nsq,png,prefect,single,svg,txt,xls
depth=1
n=1,1,1,1,2,7,1,2,3,3,1,1,1,1,2,1,1,2,2,1

>annot
annot=database_annot,go,cds,sign,bwt,soap

>exp_cov
exp_cov=coverage,rpkm(fpkm),pvalue,functional

###################
###################these are steps contained by steps above
>filter
fdir=fq
ftype=gz,gnu,png,stat,svg,xls
depth=2
n=2s,1s,2s,1s,2s,2s

>dup
fdir=fq
ftype=fq,list,stat
depth=3
n=2s,1s,2s

>as
fdir=assembly_fill
ftype=Arc,ContigIndex,R,contig,dat,edge,fa,gnu,kmerFreq,links,list,newContigIndex,peGrads,png,preArc,preGraphBasic,readInGap,readOnContig,scaf,scafSeq,scaf_gap,svg,txt,vertex,xls
step=assembly
depth=3
xdepth=2
n=1s,1s,2s,1s,3s,2s,2s,3s,1s,1s,2s,1s,1s,3s,1s,1s,1s,1s,1s,1s,1s,3s,3s,1s,2s

>f1
fdir=assembly_fill
ftype=R,dat,fa,fill,gnu,png,svg,txt,xls
step=fillgap
depth=3
xdepth=2
n=1s,2s,1s,1s,2s,2s,2s,2s,1s

>sam_clu
fdir=assembly_fill
ftype=R,Z,cidx,cluster,dat,fa,fa_cl_clusters,gnu,list,log,nhr,nin,nsq,png,singletons,svg,txt,xls
step=cluster
depth=3
xdepth=2
n=1s,1s,1s,1s,2s,5s,1s,2s,3s,3s,1s,1s,1s,2s,1s,2s,2s,1s

>database_annot
fdir=annot
ftype=Nr,R,Rout,Swissprot,cog,cut,fa_map,htm,kegg,ko,m8,path,pdf,png,xls
depth=4
n=1,1,1,1,1,1,1,1,1,1,1,1,1,1,6

>go
fdir=annot
ftype=annot,head,png,svg,wego,xls,xml
depth=2
n=1,1,1,1,1,1,1

>cds
fdir=CDS
ftype=R,conf,dat,fa,gnu,png,score,svg,txt,xls
depth=1
n=4,1,6,5,6,6,2,6,6,4

>sign
fdir=cluster
ftype=xls
depth=1
n=2

>bwt
fdir=soap
ftype=ann,bwt,fa,fmv,hot,lkt,pac,sa,sai,xls
depth=1
n=1,2,1,2,1,2,2,1,1,1

>soap
fdir=soap
ftype=PESoap,PESoapSingle
depth=1
n=1s,1s

>coverage
fdir=soap
ftype=svg,png,xls
depth=2
n=1s,1s,1s

>rpkm(fpkm)
fdir=diff
ftype=C,F,P,conf,xls
depth=1
n=1,1,1,3,3

>pvalue
fdir=diff
ftype=R,for_p,glist,ko,p,pdf,png,tmp,xls,svg,wego
depth=1,2
n=1v,1v,1v,1v,1v,1v,2v,1v,5v+3,1v,1v

>functional
fdir=diff
ftype=html,jar,png,txt
step=diff_go
depth=3
xdepth=1
n=3v+2,2,3v,3v
########################
#note:2010-12-10 change the files of coverage from R,list,pdf,png,xls to svg,png,xls
#		change the expected .xls num of CDS from 5 to 4
NUM
