use strict;
if (@ARGV < 2){
	print "perl $0 <regx, work_\*.sh> <signal, done>";
	exit()
}
my $regx = shift;
my $good_sign = shift;
my @files = glob($regx);

foreach my $i(@files){
	my @tmp = glob("$i.o*");
	# not run yet
	if (@tmp < 1){
		print "$i not_run\n";
		system("sel-qsub evo 2.5g 1 $i");

	}else{
		# run once or more, but all failed
		my $signal = 0;
		foreach my $log(@tmp){
			$signal += &check_status($log, $good_sign);
		}
		if ($signal == 0){
			print "$i run_failed\n";
			system("rm $i.e* $i.o*");
			system("sel-qsub evo 2.5g 1 $i");
		}
	}
}

sub check_status {
	my $file = shift;
	my $sign = shift;
	my $str = `cat $file`;
	if ($str =~ /$sign/){
		return 1;
	}else {
		return 0;
	}
}
