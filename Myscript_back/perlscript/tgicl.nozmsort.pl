   #!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;
use lib $FindBin::Bin;
use Mailer;

my $cap3opt = '-repeat_stringency 0.93 ';

my $usage = qq/ Usage: 
 tgicl <fasta_db> [-q <qualdb>] [-d <refDb>] [-c {<num_CPUs>|<PVM_nodefile>}]
     [-m <user>] [-O 'cap3_options']  [-l <min_overlap>] [-v <max_overhang>]
     [-p <pid>] [-n slicesize] [-s <maxsize>] [-a <cluster_file>] [-M] [-K] 
     [-L] [-X] [-I] [-C] [-G] [-W <pairwise_script.psx>] [-A <asm_program.psx>]
     [-P <param_file>] [-u <seq_list>] [-f <prefix_filter>] [-D]
 Options:
     -c : use the specified number of CPUs on local machine  ### 使用cpu数
         (default 1) or a list of PVM nodes in <PVM_nodefile>
  Clustering phase options:
     -d do not perform all-vs-all search, but search <fasta_db> against  ### 提供数据库，和数据库比对
        <refDb> instead; exit after the pairwise hits are generated
     -n number of sequences in a clustering search slice (default 1000) ### 一个cluster中 的数量
     -p minimum percent identity for overlaps <PID> (default 94)
     -l miminum overlap length (default 40)
     -G store gap information for all pairwise alignments
     -v maximum length of unmatched overhangs (default 30)
     -M ignore lower-case masking in <fasta_db> sequences
     -W use custom script <pairwise_script.psx> for the distributed ### 输入 分类的脚本
        pairwise searches instead of the default: tgicl_cluster.psx
     -Z only run the distributed pairwise searches and exit --  ### 只做 分类
       (no sorting of the pairwise overlaps and no clusters generated)
     -Y only run the distributed pairwise searches  ## 只做分类和sort 排序压缩
        and the sorted & compressed *_hits.Z file
     -L performs more restrictive, layout-based clustering
        instead of simple transitive closure
  General options:
     -I do not rebuild database indices ## 不重新 formatdb
     -s attempt to split clusters larger than <maxsize> based on  
        seeded clustering (only works if there are 'et|' 
        or 'np|'-prefixed entries provided in the input file)
     -O use given 'phrap_options' instead of the default ones ## 输出 组装参数
        ($cap3opt)
     -u skip the mgblast searches (assumed done) but restrict 
        further clustering analysis to only the sequences in <seq_list>
     -C (TIGR sequences only) always put in the same cluster all reads 
        from the same clone
     -t use <clone_list> file to put in the same cluster all sequence names 
        on the same line
     -a assemble clusters from file <cluster_file> ## 只做组装
       (do not perform any pairwise clustering)
     -f keep only sequence names with prefix <prefix> 
     -K skip the pairwise searches, only recreate the clusters
        by reprocessing the previously obtained overlaps
     -X do not perform assembly, only generate the cluster file #不组装，只生成 cluster文件
     -A use custom script as the slice assembly script 
        (instead of tgicl_asm.psx)
     -P pass the <param_file> as the custom parameter file 
        to the assembly program <asmprog.psx>
     /;

my @start_params=@ARGV;### 获取参数

my $usemail=0; #set to 1 after initialization
my $no_error=1;
my $exit_msg='';
my $cur_step;
die $usage."\n" if (@ARGV<1 || $ARGV[0] =~ /^\-/);

my $wrkdir=$ENV{'PWD'};## 获取路径

&addPath($FindBin::Bin, $FindBin::Bin.'/bin');## 将 $0加到 Bin

my $dbfile=shift(@ARGV);## 将 $ARGV[0] 作为 

umask 0002;#### 设置权限 775 

getopts('JKXSZCIDGYMLP:O:W:t:A:s:n:f:m:c:d:a:v:q:u:l:p:') || die "$usage\nError at getopts!\n";## 获取参数
my $debug=$Getopt::Std::opt_D;###debug 
my $dbfname=&getFName($dbfile);
my $dbfullpath=&getFullPath($dbfile);
my $dbqual=$Getopt::Std::opt_q;## 质量文件
my $dbqualpath;
my $clustprog=$Getopt::Std::opt_L?'lclust':'tclust';##选择clust 的程序
my $prefix=$Getopt::Std::opt_f;### 前缀
my $lstXclusive=$Getopt::Std::opt_u;### 跳过 mgblast 
die "Error: cannot find input file '$dbfile'\n" unless -s $dbfullpath;
if ($dbqual) {
	$dbqualpath=&getFullPath($dbqual);
	die "Error: cannot find quality file '$dbqual'\n" unless -s $dbqualpath;
}
my $asm_paramfile=$Getopt::Std::opt_P;
if ($asm_paramfile) {
	$asm_paramfile=&getFullPath($asm_paramfile);
	die "Error: cannot find quality file '$asm_paramfile'\n" unless -s $asm_paramfile;
}
my $cpus = $Getopt::Std::opt_c || 1;
my $usepvm = (-e $cpus) ? 1:0;
# my $useCondor= lc($cpus) eq 'condor';
my $psxcmd;
unless ($ENV{'USER'}) { # psx needs this environment variable (Solaris lacks it)
	$ENV{'USER'}=$ENV{'USERNAME'} || 'tgicl';
}

#if ($useCondor) {
#  $psxcmd="condorsx ";
#  }
# else { 
if ($usepvm) {
	$psxcmd="pvmsx -L -m $cpus ";
}else {
	die $usage."Invalid number of CPUs" unless ($cpus > 0 && $cpus <= 50);
	$psxcmd="psx -p $cpus ";
}
# }

my $mailuser = $Getopt::Std::opt_m;### 发送邮件
my $maxovh=$Getopt::Std::opt_v || 30;### maximum length of unmatched
$maxovh+=10 if $clustprog eq 'lclust';
my $minovl=$Getopt::Std::opt_l || 40;### miminum overlap length
my $pid=$Getopt::Std::opt_p || 94;### minimum percent identity for overlaps
$pid=94 if $pid<50;
$minovl=40 if $minovl==0;
my $no_zmsort=$Getopt::Std::opt_Z;###　不使用　zmsort 
my $onlysearch=$Getopt::Std::opt_Y;### only run the distributed pairwise searches and the sorted & compressed *_hits.Z file
$cap3opt=$Getopt::Std::opt_O if $Getopt::Std::opt_O;### cap3 的参数
my $clonefile;
if ($Getopt::Std::opt_t) {
	$clonefile=$Getopt::Std::opt_t;###  to put in the same cluster all sequence names on the same line
	die("Clone list file '$clonefile' not found!\n") unless -e $clonefile;
}

my $startdate=getDate();
my $nomasking=$Getopt::Std::opt_M ? 'M':'';### 忽略 序列中的 小写字母
my $gapinfo=$Getopt::Std::opt_G ? 'G':''; ## 存储gap 信息
my $useDb=$Getopt::Std::opt_d;### 仅和 user‘s db 比较
my $useDbname;
if ($useDb) {
	$useDbname=&getFName($useDb);
	$useDb=&getFullPath($useDb);
	die "Error: cannot locate the specified db: $useDb" unless -e $useDb;
}

my $psx_clust=$Getopt::Std::opt_W || 'tgicl_cluster.psx';### 选择 distributed pairwise searches的脚本
my $psx_asm=$Getopt::Std::opt_A || 'tgicl_asm.phrap.psx';### 选择 the slice assembly 的脚本

die "$psx_clust not available in PATH!\n"  unless ($psx_clust=&checkCmd($psx_clust));

die "$psx_asm not available in PATH!\n"  unless ($psx_asm=&checkCmd($psx_asm));


my $r;
#=- logging system initialization:
my $log_file="tgicl_$dbfname.log";
my $err_log="err_tgicl_$dbfname.log";
unlink($log_file, $err_log);### 删除一列文件，返回成功删除的文件名个数
open(OLDSTDERR, ">&STDERR");

open(STDOUT, ">$log_file") || &MErrExit("Failed to redirect STDOUT to $log_file");
open(STDERR, ">$err_log") || &MErrExit("Failed to redirect STDERR to $err_log");
&set_step('Initialization');### 输出日志，记录当前步骤
### &flog等于 print STDOUT 和 print STDERR ;
&flog("$FindBin::Script running options: \n".$FindBin::Script.' '.join(" ", @start_params));
&flog(" Standard log file: $log_file");
&flog(" Error log file:    $err_log");
&flog(" Using $cpus CPUs for clustering and assembly");
&flog(" Path is : $ENV{'PATH'} ");
### 输出日志 记录当前 使用的参数 和 当前系统变量
$no_error=0;
$usemail=1 && $mailuser;

my $clusterfile=$dbfname.($useDb? '_'.$useDbname :'_cl').'_clusters';### 一个的输出文件命名

my ($cmd, $psxparam);
$Getopt::Std::opt_I=1 if $Getopt::Std::opt_K;###　skip the pairwise searches, only recreate the clusters　by reprocessing the previously obtained overlaps　

my $filterfile=$dbfname.'_filter.lst';### 一个的输出文件命名

#
# rebuild indices (unless specifically told not to do so)
unless ($Getopt::Std::opt_I) {
	my $toindex = $useDb ? $useDb : $dbfile ;
	&flog("-= Rebuilding $toindex indices =-");
	$cmd="formatdb -p F -o F -i $toindex";
	system($cmd) &&  &MErrExit("Error at '$cmd'");### system 返回退出状态 ，正常结束 返回 0 
	system("cdbfasta $toindex") && &MErrExit("Error at cdbfasta $toindex");
	if ($dbqual && !-e $dbqual.'.cidx') {  
		system("cdbfasta $dbqual") &&   &MErrExit("Error at cdbfasta $dbqual");   
	}
}

my $dbfileidx=$dbfullpath.'.cidx';

if ($Getopt::Std::opt_S) {### 这个太无赖了，-S 就没在 usage 中出现。。
	$clusterfile=$Getopt::Std::opt_a if $Getopt::Std::opt_a;####直接只做组装　通过　-a 输入
	goto SINGLETSONLY;
}

if ($Getopt::Std::opt_a) {####直接只做组装　通过　-a 输入
	$clusterfile=$Getopt::Std::opt_a;
	&flog("-- Skipping clustering phase.\n");
	goto ASSEMBLE;
}

CLUSTER:
#=- start clustering the $blastDB file 开始聚类分析
&set_step('clustering');
unlink($clusterfile);
if ($Getopt::Std::opt_K) {##skip the pairwise searches, only recreate the clusters by reprocessing the previously obtained overlaps
	&flog("-- Skipping pairwise searches.\n");
	goto TCLUSTER;
}
my $numseqs = $Getopt::Std::opt_n || 1000;### number of sequences in a clustering search slice 
system('/bin/rm -rf cluster_[1-9]*') unless $Getopt::Std::opt_J;### 这个太无赖了，-J 也没在 usage 中出现。。
# ---psx user parameter format:   <db>:<minpid>:<maxovh>:<minovl>:<flags>
#  where <flags> are one or more letters of: 
#      D=no self-clustering, M = no masking, G = gap info
my $dbflag=$useDb?'D':'';
my $paramdb=$useDb? $useDb : $dbfullpath;
$psxparam=join(':', ($paramdb,$pid,$maxovh,$minovl,$nomasking.$gapinfo.$dbflag));
$cmd=$psxcmd." -n $numseqs -i $dbfile -d cluster -C '$psxparam' -c '$psx_clust'";
unless ($Getopt::Std::opt_J) {
	&flog("  Launching distributed clustering: \n $cmd");
	system($cmd)  && &MErrExit("Error at '$cmd'\n");
}else { &flog(" Skipping the mgblast searches.. (-J option given)"); }
&end_step();

#========= compress & merge sort the cluster results

if ($no_zmsort) {
	&flog("Exit requested after pairwise searches.");
	# &flog($exit_msg);
	# generate the clustering singleton list here?
	goto THEEND;
}



my @dirs = ((<cluster_?>),(<cluster_??>));### 获取cluster 的路径 cluster_? 和 cluster_?? 
unlink('masked.lst');##删除 'masked.lst'
foreach my $dir (@dirs) {
	next unless -d $dir;
	system("cat $dir/masked.lst >> masked.lst");
	#$cmd="zmsort -f11 -n -r -o zdir_$dir -s 700 $dir/*.Z";
	zmergeDirHits($dir);
	#system($cmd)
	#   && &MErrExit("Error at command:\n$cmd");
	#system("/bin/rm -rf $dir") unless $debug;
}
my $hitsort = $dbfname.'_'.($useDb ? $useDbname.'_tabhits' : 'cl_tabhits');
#$cmd="zmsort -f10 -n -r -o $hitsort -s 1600 zdir_*.Z";### 将每个 cpu 产生 zmsort 文件 合并
$cmd="gzip -cd  zdir_*.Z |sort -k 10 -n -r -o $hitsort && gzip $hitsort && mv $hitsort.gz $hitsort\_1.Z";### add by nixiaoming
system($cmd)   && &MErrExit("Error at final sort:\n $cmd");
system('/bin/rm -rf zdir_*.Z');

if ($onlysearch) {### 如果仅作 sort ，那么到此处 终止
	&flog("Exit requested after pairwise searches and sorting");
	# &flog($exit_msg);
	# generate the clustering singleton list here?
	goto THEEND;
}

TCLUSTER:
my %lst;
if ($lstXclusive) {### 是否跳过 mgblast 
	open(SEQLST, $lstXclusive) || die "Error opening sequence list '$lstXclusive'\n";
	while (<SEQLST>) {
		next if m/^\s*#/;
		chomp;
		next unless $_;
		$lst{$_}=1;
	}
	close(SEQLST);
}

if ($Getopt::Std::opt_C) { #prepare clone list if -C ## (TIGR sequences only) always put in the same cluster all reads from the same clone
	if ($prefix || $lstXclusive) {### -f 和 -u 
		open(SEQONLY, '>'.$filterfile) || die "Error creating file '$filterfile'";### 这句似乎有错 应该是 '>', 而不是 . ## $filterfile=$dbfname.'_filter.lst';### 一个的输出文件命名
	}
	$cmd = "cdbyank -l $dbfileidx | sort |"; ### $dbfullpath.'.cidx'; ## 这个将所有的 序列名称提取出来了 
	$clonefile = $dbfname.'.clones';## 第128行  $clonefile=$Getopt::Std::opt_t;###  to put in the same cluster all sequence names on the same line

	open(CLONES, '>'.$clonefile) || die "Error creating file '$clonefile'";## 新发现，这个尽可以 使用 '>'. ## 这个文件似乎是 按序列名称的 前 7个 字符 分类，这不适用 与 我们转录组的情况 
	#&flog("Fetching output of '$cmd'");
	open(RLIST, $cmd);####还可以这么用？在open中 使用命令行 ，注意要以 管道 | 结尾
	my @clonereads;
	while (<RLIST>) {
		chomp;
		s/^\s+//;
		if ($prefix) {
			next unless $prefix eq substr($_,0,length($prefix));
			print SEQONLY $_."\n";
		}elsif ($lstXclusive) {
			next unless exists $lst{$_};
			print SEQONLY $_."\n";
		}
		my $clone=substr($_,0,7);###这个　7 是怎么定的？ 是否不适用于 我们转录组的情况 unigene scaffold 这个都超过了 7个字符..
		if ((@clonereads==0) || $clone eq substr($clonereads[0],0,7)) {
			push(@clonereads, $_);
		}else  {
			print CLONES join(' ',@clonereads),"\n" if @clonereads>1;
			@clonereads=($_);
		}
	}
	print CLONES join(' ',@clonereads),"\n" if @clonereads>1;
	close(RLIST);
	close(CLONES);
	$cmd="gzip -cd ${hitsort}_*.Z | $clustprog PID=$pid OVL=$minovl OVHANG=$maxovh -o $clusterfile -c $clonefile";### my $clustprog=$Getopt::Std::opt_L?'lclust':'tclust';##选择clust 的程序 ## my $pid=$Getopt::Std::opt_p || 94;### minimum percent identity for overlaps ## my $minovl=$Getopt::Std::opt_l || 40;### miminum overlap length ## my $maxovh=$Getopt::Std::opt_v || 30;### maximum length of unmatched

	if ($prefix || $lstXclusive) {
		close(SEQONLY);
		$cmd.=" -r $filterfile";
	}
}else { # do not generate TIGR clone list  
	$cmd="gzip -cd ${hitsort}_*.Z | $clustprog PID=$pid OVL=$minovl OVHANG=$maxovh -o $clusterfile ";
	$cmd.= "-c $clonefile"  if $clonefile;
	if ($prefix) {
		my $pcmd="cdbyank -l $dbfileidx | grep '^$prefix' > $filterfile";
		system($pcmd) && &MErrExit("Error building filtered sequence name list:\n$pcmd\n");
		$cmd.=" -r $filterfile";
	}elsif ($lstXclusive) {
		my $pcmd="cat $lstXclusive > $filterfile";
		system($pcmd) && &MErrExit("Error building filtered sequence name list:\n$pcmd\n");
		$cmd.=" -r $filterfile";
	}
}

&flog("Running transitive closure command: $cmd\n");
system($cmd)  && &MErrExit("Error running transitive closure step:\n$cmd\n");
if ($Getopt::Std::opt_s) { #try to break too large clusters based on full mRNA 'et|' prefixed entries
	## 转录组的可以考虑 将长度大于一定长度的作为 全长转录本 比如 5000 
	my $clmax=$Getopt::Std::opt_s;### 将长度大于一定长度 ( > $clmax )的作为 全长转录本
	&flog("Attempting to split clusters with more than $clmax sequences");
	open(CLFILE, $clusterfile) || &errdie("Cannot open cluster file $clusterfile\n");## 打开聚类文件
	open(LARGECLS, '>'.$clusterfile.'.large') || &errdie("Cannot create $clusterfile.large\n");
	open(SMALLCLS, '>'.$clusterfile.'.small') || &errdie("Cannot create $clusterfile.small\n");
	my $islarge=0;
	while (<CLFILE>) {
		if (m/^>\w+\s+(\d+)/) {
			$islarge = ($1 > $clmax);
		}
		if ($islarge) {
			print LARGECLS $_;
		}else {
			print SMALLCLS $_;
		}
	}
	close(LARGECLS);close(SMALLCLS); close(CLFILE);
	if (-s $clusterfile.'.large') {### 如果存在 .large 将 .large 从 $clusterfile 中剔除，并对 .small 进行 sclust 处理
		#create the xclusion file for sclust:
		system("grep -v '^>' $clusterfile.small | tr ' \\t' '\\012\\012' > xclude_small") && &errdie("Error collecting small clusters sequence list");
		system("gzip -cd ${hitsort}_*.Z | sclust -L -x xclude_small HEAVY=3200 -o largecls.scls")  && &errdie("Error running sclust!");
		system("cat largecls.scls $clusterfile.small > $clusterfile")  && &errdie("Error at 'cat largecls.scls $clusterfile.small > $clusterfile'");
		unlink('largecls.scls', 'xclude_small');
	}
	unlink($clusterfile.'.small', $clusterfile.'.large');
}
&flog("The clusters are stored in file '$clusterfile'.\n");

if ($Getopt::Std::opt_X || $useDb) {#### -X 只做聚类，不做组装
	&flog("Exit requested after clustering phase.");
	# &flog($exit_msg);
	# generate the clustering singleton list here?
	goto THEEND;
}
  
ASSEMBLE:
&set_step('ASSEMBLE');
&MErrExit("Error: cluster file '$clusterfile' is missing or empty.") unless (-s $clusterfile);
my $dbqualidx;
if ($dbqual) {
	$dbqualidx=$dbqualpath.'.cidx';
	&MErrExit("Error: cdbfasta index ($dbqualidx) cannot be located.")  unless (-e $dbqualidx);
}
 
&MErrExit("Error: cdbfasta index ($dbfileidx) cannot be located.") unless (-e $dbfileidx);

system("/bin/rm -rf asm_[1-9]*");
$cap3opt =~ s/\s+/\~\+/g;
$psxparam = $dbfileidx.':'.$cap3opt.':'.$dbqualidx;
$psxparam.=':'.$asm_paramfile; 
$psxparam.=":$pid~$minovl~$maxovh";
$cmd="$psxcmd -i $clusterfile -d asm -C '$psxparam' -c '$psx_asm'"; ##$psxcmd=psx -p $cpus
system($cmd);
if ($?) {
	$exit_msg=`cat asm_*/err_log`;
	&MErrExit("Error at: \n$cmd\n");
}
my $errlog=`cat asm_*/err_log`;
$errlog =~ s/^\s+//g;
($errlog =~ s/An error occurred in construction of a consensus//g);
$errlog =~ s/ALIGNMENT ERROR/ALIGNMENT/g;### add bu nixiaoming # phrap : ALIGNMENT ERROR Unigene49279_FBH_ONCfkhTBRAAPEI-9 -174 -5  509 Done

if ($errlog =~ m/error/is) {
	$exit_msg="Errors detected in the assembly subdirectories:\n".$errlog;
	exit(0);
}

SINGLETSONLY:
unless ($Getopt::Std::opt_a) { ### -a 只做组装
	#build the singleton list if the full straight clustering/assembly pipeline was run
	$cmd='cat asm_*/ACE | grep "^AF " | cut -d " " -f2 | sort -u > seqnames_in_asms';
	system($cmd) && &MErrExit("Error running:\n$cmd\n");
	if (($prefix || $lstXclusive) && -s $filterfile) {
		$cmd="sort -o seqnames_all $filterfile";
	}else {
		$cmd="cdbyank -l $dbfileidx | sort > seqnames_all";
	}
	system($cmd) && &MErrExit("Error running:\n$cmd\n");
	my $sglist = $dbfname.'.singletons';
	$cmd="comm -23 seqnames_all seqnames_in_asms > $sglist";### 获取未被聚类的序列id
	system($cmd) && &MErrExit("Error running:\n$cmd");
	unlink('seqnames_all', 'seqnames_in_asms');
}
  
&end_step();

THEEND:
$no_error=1;
&flog("*** tgicl [$dbfile] finished ***");


#=====================================================================
#=========================    SUBROUTINES   ==========================
#=====================================================================

END { #to be executed on exit
	if ($cur_step && $err_log) {
		my $step_log=`cat $err_log`;
		$step_log =~ s/\s+$//;
		$exit_msg.="\nThis is the content of the error log file (ended at $cur_step):\n$step_log"
		if $step_log;
		my $host=$ENV{'HOST'} || $ENV{'HOSTNAME'};
		my $msg = $no_error ? qq/$FindBin::Script ($dbfile) finished on machine $host
                 in $wrkdir, without a detectable error.
                 / :
                 qq/$FindBin::Script ($dbfile) encountered an error at step $cur_step
                 Working directory was $wrkdir.                 
                 /;
		unless ($no_error) { #an abnormal termination
			&flog("\nProcess terminated with an error, at step '$cur_step'!");
			&send_mail({to=>$mailuser, subj=>"$FindBin::Script ($dbfile) error at $cur_step!",body=>$msg.$exit_msg}) if $usemail;   
			&flog($msg);
		}else {
			#&flog("*** Done ***") if ($cur_step && lc($cur_step) ne 'Initialization');
			&send_mail({to=>$mailuser, subj=>"$FindBin::Script ($dbfile) finished.",body=>$msg.$exit_msg}) if $usemail;
		}
		print OLDSTDERR $msg;
	}
}

#== checkCmd -- checks for executable, in the PATH if no full path given
sub checkCmd {
	my $cmd=$_[0];
	if ($cmd =~ m/^\//) {
		return (-x $cmd) ? $cmd : '';
	}
	my @paths=split(/:/, $ENV{'PATH'});
	foreach my $p (@paths) {
		return $p.'/'.$cmd if -x $p.'/'.$cmd;
	}
	return ''; 
}

sub flog {
	print STDOUT join("\n",@_),"\n";
	print STDERR join("\n",@_),"\n";
}

sub MErrExit {
	#print STDERR $_[0]."\n";
	$exit_msg.=$_[0].$_[1];
	&flog($exit_msg);
	exit(1) unless defined($_[1]);
	die $_[1];
}


sub set_step {
	$cur_step=$_[0];
	&flog(">>> --- $cur_step [$dbfile] started at ".&getDate());
}

sub end_step {
	&flog("<<< --- $cur_step [$dbfile] finished at ".getDate());
}

#a basic date function :
sub getDate {
	my $date=localtime();
	#get rid of the day so Sybase will accept it
	(my $wday,$date)=split(/\s+/,$date,2);
	return $date;
}

sub getFullPath {
	return ($_[0] =~ m/^\//) ? $_[0] : $ENV{'PWD'}.'/'.$_[0];
}

sub getFName {
	if ($_[0] =~ m/.+[\/\\](.*?)$/) {
		return $1;
	}else {
		return $_[0];
	}
}

sub addPath {
	my $path=$ENV{'PATH'};
	foreach my $p (@_) {
		next if ($p eq $path || m/\Q:$p$/  || m/:\Q$p:/ || m/^\Q$p:/);
		$path=$p.':'.$path;
	}
	$ENV{'PATH'}=$path; 
}

#sub zmergeDirHits {
#	my ($dir)=@_;
#	my $maxfopen=20; #max 20 files at once ### 可以考虑 改小
#	my $run=0;
#	chdir($dir) || die "Cannot change to $dir directory!\n";
#	while (1) {
#		opendir(FDIR, '.') || die "Cannot open directory $dir\n";
#		my @allfiles=readdir(FDIR); #possibly large array w/ all the files from that directory 将$dir 下的所有文件名读入到 @allfiles , 包括 . .. .* 
#		close(FDIR);
#		@allfiles=grep(/\.Z$/, @allfiles);### 提取 所有 .Z 的 文件
#		last if (@allfiles<=$maxfopen);### 如果 文件总数 过少则 跳出
#		my @files;
#		foreach my $f (@allfiles) {
#			next unless ($f=~m/\.tab\.Z$/ || $f=~m/zMrg_p\S+\.Z$/);### zMrg_p*.Z 是输出文件 .tab.Z 是 所要的输入文件 
#			push(@files, $f);
#			if (@files==$maxfopen) {
#				#my $sortcmd="mgmerge -o zMrg_p$run -s 1200 ".join(' ',@files);
#				my $sortcmd="zmsort -f10 -n -r -o zMrg_p$run -s 1200 ".join(' ',@files);
#				$run++;
#				&runZmerge($sortcmd);### 执行 $sortcmd 并检查 返回状态，输出日志 到 STDERR
#				if ($debug) {
#					foreach my $rf (@files) {
#						my $ren=$rf;
#						$ren=~s/\.Z$/\.gz/;
#						rename($rf, $ren); ### rename 函数 是对文件重命名 相当于 mv 
#					}
#				}else {
#					unlink(@files);### unlink 函数 是删除 文件，相当于 rm 
#				}
#				@files=();## 清空 @files 
#			}
#		}
#		if (@files>0) { ### 文件数少于 20 的部分进行 处理，是上个循环的扫尾 
#			#my $sortcmd="mgmerge -o zMrg_p$run -s 1200 ".join(' ',@files);
#			my $sortcmd="zmsort -f10 -n -r -o zMrg_p$run -s 1200 ".join(' ',@files);
#			$run++;
#			&runZmerge($sortcmd);
#			if ($debug) {
#				foreach my $rf (@files) {
#					my $ren=$rf;
#					$ren=~s/Z$/gz/;
#					rename($rf, $ren);
#				}
#			}else { 
#				unlink(@files); 
#			}
#		}
#	};
#	chdir('..');
#	my $sortcmd="zmsort -f10 -n -r -o zdir_$dir -s 1600 $dir/*.Z";### zmsort 处理 while 循环产生 的 zMrg_p$run.Z 
#	runZmerge($sortcmd);
#	system("/bin/rm -rf $dir") unless $debug;
#}

sub zmergeDirHits {#### 因为zmsort 在文件行数 超过10000 时报错，改用 sort 试试
	my ($dir)=@_;
	chdir($dir) || die "Cannot change to $dir directory!\n";
	opendir(FDIR, '.') || die "Cannot open directory $dir\n";
	my @allfiles=readdir(FDIR); #possibly large array w/ all the files from that directory 将$dir 下的所有文件名读入到 @allfiles , 包括 . .. .* 
	close(FDIR);
	@allfiles=grep(/\.tab\.Z$/, @allfiles);### 提取 所有 .Z 的 文件
	my $tmp2='';
	for (my $i=0;$i*20 <scalar @allfiles ;$i++) {
		my $tmp='';
		for (my $j=0;$j<20 ;$j++) {
			if ($i*20+$j<scalar @allfiles) {
				$tmp.=' '.$allfiles[$i*20+$j];
			}else{
				last;
			}
		}
		my $sortcwd="/bin/gzip -cd $tmp |sort -n -r -k 10 -o sort.tmp.$i.xls ";
		$tmp2.=" sort.tmp.$i.xls";
		system "$sortcwd" ;#&& &MErrExit("Error at  gzip:\n $sortcwd ");
	}
	my $sortcwd="cat  $tmp2 |sort -n -r -k 10 -o  ../zdir_$dir  ";
	system "$sortcwd" ;#&& &MErrExit("Error at  gzip:\n $sortcwd ");
	if ($debug) {
		foreach my $rf (@allfiles) {
			my $ren=$rf;
			$ren=~s/Z$/gz/;
			rename($rf, $ren);
		}
	}else { 
		unlink(@allfiles); 
	}
	chdir('..');
	#system "sort -k 10 -n -r $dir/sort.tmp.xls -o zdir_$dir " && &MErrExit("Error at  sort:\n sort $dir/sort.tmp.xls ");
	system "gzip zdir_$dir ";# && &MErrExit("Error at  sort:\n gzip zdir_$dir ");
	rename("zdir_$dir.gz","zdir_$dir.Z");
	system("/bin/rm -rf $dir") unless $debug;
}
#sub runZmerge {
#	my $cmd=shift(@_);
#	print STDERR "Running:\n$cmd\n" if $debug;
#	my $errout=`( $cmd ) 2>&1`;
#	my $errcode=$?;#### 这个值记录了上一次管道关闭时的返回状态，`` 或 system 或 wait 或waitpid； 成功时 返回 0
#	print STDERR $errout;
#	if ($errcode) {
#		print STDERR "Error status detected (code=$errcode) at command $cmd\n";
#		exit(1); #exit, but do not report error 
#		 #- so we can continue to assemble the rest of the clusters!
#	}
#	if ($errout=~/aborted|error|fail|segmentation|violation|cannot/i) {
#		print STDERR "Error message detected after running command:\n$cmd\n";
#		exit(2);
#	}
#}
