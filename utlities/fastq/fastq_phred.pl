#!/usr/bin/perl -w
use strict;
use Getopt::Long;
#
# fastq_phred.pl - Script for judge the fastq's encoding, whether it is phred33 or phred64.
#
# Version: 0.3 ( May 19, 2014)
# Author: Wencai Jie (jiewencai<@>qq.com), NJAU, China.
#
# Permission is granted to anyone to use this software for any purpose, without
# any express or implied warranty. In no event will the authors be held liable 
# for any damages arising from the use of this software.
#

#Get options.
my ($help, $print_score, $detail, $print_ascii, $reads_num, $reads_start_arg, $reads_end_arg);
my $reads_end_turn;
GetOptions(
	'help|h!' => \$help,
	'score|s!' => \$print_score,
	'detail|d!' => \$detail,
	'ascii|a!' => \$print_ascii,
	'reads_num|n=i' => \$reads_num,
	'reads_start|b=i' => \$reads_start_arg,
	'reads_end|e=i' => \$reads_end_arg,
);

my $usage = "
fastq_phred.pl:
This program can print fastq file's reads quality scores, ASCII value, and help to judge it's 
encoding by the ASCII value range, whether it is phred33 or phred64.

Usage:
	perl fastq_phred.pl [options] <file1.fq [file2.fq ...]>
Options:
	-h|--help         print this help message.
	-s|--score        print scores.                                 [default: Do not print scores] 
	-d|--detail       print detail scores or ASCII value when       [default: Do not print scores in detail] 
                          --score or --ascii set.
	-a|--ascii        print quality character's ASCII value. if     [default: Do not print ASCII vaule] 
                          this option set, the --score will disabled.
	-n|--reads_num    reads number used to test phred encoding      [default: 1000]
                          and print scores. It's advised to use more
                          than 100 reads to do the test.               
	-b|--reads_start  reads start position used to test phred       [default: 1]
                          encoding and print scores.
	-e|--reads_end    reads end position used to test phred         [default: the length of the read]
                          encoding and print scores.

";
	
if ($#ARGV < 0 or $help){ 
	print "$usage";
	exit;
}

#Check parameters.
unless ($reads_num){
	$reads_num = 1000;
}
if ($reads_start_arg && $reads_start_arg  < 0){
	print STDERR "ERROR:The reads start position should great than 0.\n\n";
	exit;
}
if ($reads_end_arg && $reads_end_arg  < 0){
	print STDERR "ERROR:The reads end position should great than 0.\n\n";
	exit;
}
if ($reads_start_arg && $reads_end_arg && $reads_end_arg < $reads_start_arg){
	print STDERR "ERROR:The reads start position should great than end position.\n\n";
	exit;
}

&main;

sub main{
	my $filename = '';
	while ($filename = shift @ARGV){
	my @FQ = ();
	my @all_ascii = ();
	my ($file_end, $phred_result) = ('','');
	my ($Q, $count, $lt_58, $gt_75) = (0, 0, 0, 0);
    if ($filename =~ /.gz$/) {
	    open FQ,"<:gzip","$filename" or die "Can not open $filename:$!\n";
    }else{
        open FQ,"$filename" or die "Can not open $filename:$!\n";
    }
	#Read sequences.
	while($count < $reads_num){
		$count++;
		@FQ=();	
		#read four lines from fastq file.
		for(my $i=0; $i<=3; $i++){
			if (eof(FQ)){
				$file_end = 'yes';
				last;
			}
			$FQ[$i]=<FQ>;
			if ($FQ[0] !~ m/^@/){
				my $line = $count*4-3;
				print STDERR "ERROR:\n$filename: It's not a correct fastq format.\nline '$line': $FQ[0]\n";
				exit;
			}
		}
		if ( $file_end eq 'yes'){
			next;
		}
		my @ascii_ref = &cal_ascii($FQ[3], $reads_start_arg, $reads_end_arg);
		push @all_ascii, [@ascii_ref];
	}

	#print ASCII.
	if ($print_ascii){
		print "\n","."x50," ASCII Value: $filename ","."x50,"\n";
		&print_array_of_array(\@all_ascii, 0, $detail);
		next;
	}

	#Stastic ASCII value range.
	foreach my $ascii_ref (@all_ascii){
		$lt_58 += (grep { $_ <= 58} @{$ascii_ref});
		$gt_75 += (grep { $_ >= 75} @{$ascii_ref});
	}

	#Guess the Phred with ASCII value range. 
	if ($lt_58 > 1 && $gt_75 == 0 ){
		$Q = 33;
		$phred_result = "$filename: The encoding should be Phred33.\nThe quality score character number that ASCII value less than 58 : $lt_58\nThe quality score character number that ASCII value great than 75: $gt_75";
	}elsif($lt_58 == 0 && $gt_75 > 1){
		$Q = 64;
		$phred_result = "$filename: The encoding should be Phred64.\nThe quality score character number that ASCII value less than 58 : $lt_58\nThe quality score character number that ASCII value great than 75: $gt_75";
	}elsif($lt_58 == 0 && $gt_75 == 0){
		print STDERR "$filename: The encoding should be Phred33 that all of the nucleotide quality score great than 25 and less than 41, but it's advised to send more reads to be tested with '-n <int>' options.\n";
		exit;
	}else{
		print STDERR "$filename\nWarning: Abnormal endoding, Please test again with more reads or make a judgement by yourself with ASCII value by '-ascii' options.\n"; 
		exit;
	}

	#print score.
	if ($print_score){
		print "\n","."x50," Quality Score: $filename ","."x50,"\n";
		&print_array_of_array(\@all_ascii, $Q, $detail);
	}

	#Print the phred encoding result.
	print STDERR "$phred_result\n\n";
	}
}

#Print Score or ASCII value.
sub print_array_of_array{
	my ($array_of_array_ref, $Q, $detail) = @_;
	my ($average_value, $total_value, $value_num) = (0, 0, 0);
	my @array_of_array = @{$array_of_array_ref};
	my %value_h;
	foreach my $array_ref (@array_of_array){
		for (my $i=0;$i<=$#{$array_ref};$i++){
			my $out_value  = ${$array_ref}[$i] - $Q;
			$value_h{$out_value}++;
			$total_value += $out_value;
			$value_num ++;
			print "$out_value " if ($detail);
		}
		print "\n" if ($detail);
	}
	unless ($detail){
		foreach my $out_value (sort {$a <=> $b} keys %value_h){
			print "$out_value\t$value_h{$out_value}\n";
		}
	}
	$average_value = (int ($total_value/$value_num)*100) /100;
	print "Average: $average_value\n";
}

#Calculate phred score.
sub cal_ascii{
	my ($read,$reads_start, $reads_end) = @_;
	my @all_ascii = ();
	my $ascii = 0;
	#The $read string's end is a "\n";
	my $reads_len = length($read) - 1;
	#The $reads_end should be less than the read length.
	if( $reads_end_arg && $reads_end_arg <= $reads_len){
		$reads_end = $reads_end_arg;
	}else{
		$reads_end = $reads_len;
	}
	#If the the reads start position set, the $reads_start equal to it,
	#else the reads start position set to 1.
	if ( $reads_start_arg && $reads_start_arg <= $reads_end){
		$reads_start = $reads_start_arg;
	}else{
		$reads_start = 1;
	}
	#Convert 1 base coordinate system to 0 base coordinate system.
	for(my $j=$reads_start-1; $j<=$reads_end-1; $j++){
		$ascii = ord(substr($read,$j,1));
		push @all_ascii, $ascii;
	}
	return @all_ascii;
}
