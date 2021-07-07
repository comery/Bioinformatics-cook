#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage=<<_USAGE_;
Description

    to monitor qsub task and then kill it or not.

Usage

    perl $0  -j <job_ID>  [-t over_time]  [-i check_interval]  [-k]  [-h]

    -j <INT>   qsub job ID.
    -t [INT]   the lasting time (minutes) since the job exceed the #virtual_free, i.e. the #vf value when qsub the job. default: 5.
    -i [INT]   the frequency of checking (minutes). default: 5.
    -k         kill the job if it exceeds the #virtual_free for -t <INT> minutes.
    -h         this infor.

_USAGE_


my ($job, $kill, $help);

my $overTime = 5;
my $interval = 5;

GetOptions(
    "j:i"  =>  \$job,
    "t:i"  =>  \$overTime,
    "i:i"  =>  \$interval,
    "k!"   =>  \$kill,
    "h|help!" => \$help,
);

die "$usage\n" unless (defined $job);

$overTime *= 60;
$interval *= 60;

my $msg = `qstat -j $job | sed -n "17p"`;
$msg =~ /^hard\s+resource_list\:\s+virtual_free\=(.+)G$/i || die "failed to get #virtual_free!\n";
my $virtual_free = $1;

sleep(60);

my ($i, $time, $vmem);
open LOG, ">$job.qstat.log";
my $old_fh = select(LOG);
$| = 1;

while (1){
	$time = localtime();
	$msg = `qstat -j $job | sed -n "26p"`;
    
    print "$time\nvirtual_free=${virtual_free}G\n$msg";

    if($msg=~/^usage.+vmem\=(.+?[GM]), maxvmem\=(.+[GM])$/){
	    $vmem = $1;
        if ($vmem=~s/M$//){
            $vmem = $vmem/1000;
        }elsif ($vmem=~s/G$//){
            $vmem = $vmem;
        }
    }elsif($msg=~/queue\s+instance/){ ## not begin running yet
        print "This job is waiting to start...\n";
        sleep($interval);
        next;
    }else{
        print "This job maybe has been finised!\n";
		last;
	}
	
	if ($kill){
		&kill_job();
	}

	sleep($interval);
}
	
select($old_fh);
    
sub kill_job{
	if ($vmem > $virtual_free){ 
		print "I am going to kill this job if the #vmem > ${virtual_free}G for $overTime seconds...\n";
		$i++;

		if (($interval*$i) > $overTime){
		    `qdel $job`;
			print "The job $job has been killed!\n\n";
		    last;
		}

	}else{
		$i = 0;
    }
}

