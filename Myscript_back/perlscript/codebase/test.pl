#!/usr/bin/perl -w 

=head1 Name
	
	split.pl -split high-throughout barcoding data into different file with index

=head1 Usage
	
	perl $0  [option]
	--1 <str> fastq1
	--2 <str> fastq2
	--primer <str> index_table.xls
	--outdir <str> output dirname
	--l <number> the length of index+primer
	--type <str> fa|fq|both

=head1 Exmple

	perl split_by_index.pl -1 hifi_1.fq -2 hifi_2.fq -primer index_table.xls -outdir index -l 30

=cut
use strict;
use Getopt::Long;
my ($left,$right,$primer,$outdir,$cut_len,$help);
GetOptions (
		"1:s" =>\$left,
		"2:s" =>\$right,
		"primer:s" =>\$primer,
		"outdir:s" =>\$outdir,
		"l:i" =>\$cut_len,
		"help" =>\$help
);
die `pod2text $0` if (!($left && $right && $primer && $cut_len)|| $help);
$outdir ||= "index";


#-----------remove old  outdir-------
if (! -d $outdir) {
	mkdir ("$outdir") 
}else{
	`rm -r $outdir/*`;
}
#-----------------------------------
open INDEX,$primer;
my (%index1,%index2,%index3,%index4);
while (<INDEX>) {
	chomp;
	next if ($_ =~ /^num/);
	my @a = split /\s+/,$_;
	my $num = $a[0];
	`mkdir -p $outdir/$num`;
	my $for = $a[1];
	my $rev = $a[2];
	my $revFor = &revcom($for);
	my $revRev = &revcom($rev);
	$index1{$for} = "For".$num;
	$index1{$rev} = "Rev".$num;
	$index2{$revFor} = "revFor".$num;
	$index2{$revRev} = "revRev".$num;
}
close INDEX;

open LE, $left;
open RI, $right;
open MID1, ">middle_1.fa";
open MID2, ">middle_2.fa";
open MID3, ">middle_1.fq";
open MID4, ">middle_2.fq";
open ERR1, ">err_1.fq";
open ERR2, ">err_2.fq";
my $count_op = 0;
my $count_ng = 0;
#-----------read fastq file---------------
$/="@";<LE>;<RI>;$/="\n";
while (my $id_le = <LE>,my $id_ri = <RI>) {
	chomp ($id_le,$id_ri);
	$/="@";
	my $stuff_le = <LE>;
	my $stuff_ri = <RI>;
	chomp $stuff_le;chomp $stuff_ri;
	my $seq_le = &seq($stuff_le);
	my $seq_ri = &seq($stuff_ri);
	print "fastq1 and fastq2 are not in order!" if ($id_le ne &pair_name($id_ri)) ;
	my $head_1 = substr ($seq_le,0,$cut_len);
	my $head_2 = substr ($seq_ri,0,$cut_len);
	my $tail_1 = substr ($seq_le,-$cut_len);
	my $tail_2 = substr ($seq_ri,-$cut_len);
	my $file;
	if (exists $index1{$head_1} && ! exists $index2{$head_2} || exists $index1{$head_2} &&  ! exists $index1{$head_2})  {     # positive direction "+"
		$count_op ++;
		$file = &output($index1{$head_1});
	#	print "$file\n";
		open OUT1, ">>$file\_1.fa";open OUT2, ">>$file\_2.fa";
		open OUT3, ">>$file\_1.fq";open OUT4, ">>$file\_2.fq";
		print OUT1 ">"."$id_le\n$seq_le\n";
		print OUT2 ">"."$id_ri\n$seq_ri\n";
		print OUT3 "@"."$id_le\n$stuff_le";
		print OUT4 "@"."$id_ri\n$stuff_ri";
		close OUT1;close OUT2;close OUT3,close OUT4;
	}elsif (exists $index1{$head_1} && exists $index2{$head_2} || exists $index1{$head_2} && exists $index1{$head_2})  { 
		print ERR1 "@"."$id_le\n$stuff_le";  #error
		print ERR2 "@"."$id_ri\n$stuff_ri";  #error

	}elsif (exists $index2{$head_1} && ! exists $index2{$head_2} || exists $index2{$head_2} && ! exists $index2{$head_1}) {  # negative direction "-"
		$count_ng ++;
		$file = &output($index1{$head_2});
		open OUT1, ">>$file\_1.fa";open OUT2, ">>$file\_2.fa";
		open OUT3, ">>$file\_1.fq";open OUT4, ">>$file\_2.fq";
		print OUT1 ">"."$id_le\n".&revcom($seq_le)."\n";
		print OUT2 ">"."$id_ri\n".&revcom($seq_ri)."\n";
		print OUT3 "@"."$id_le\n".&revfq($stuff_le)."\n";
		print OUT4 "@"."$id_ri\n".&revfq($stuff_ri)."\n";
		close OUT1;close OUT2;close OUT3,close OUT4;
	}elsif (exists $index2{$head_1} && exists $index2{$head_2} || exists $index2{$head_2} && exists $index2{$head_1}) { 
		print ERR1 "@"."$id_le\n$stuff_le";  #error
		print ERR2 "@"."$id_ri\n$stuff_ri";  #error
	
	}elsif (exists $index2{$tail_1}  || exists $index2{$tail_2}) { # primer in middle of reads
		print ERR1 "@"."$id_le\n$stuff_le";  #error
		print ERR2 "@"."$id_ri\n$stuff_ri";  #error
	
	}else {
		print MID1 ">"."$id_le\n$seq_le\n";
		print MID2 ">"."$id_ri\n$seq_ri\n";
		print MID3 ">"."$id_le\n$stuff_le";
		print MID4 ">"."$id_ri\n$stuff_ri";
	}

	$/= "\n";
}
print "$count_op\t$count_ng";
close MID1;
close MID2;
close ERR1;
close ERR2;
#===================================================
#               SUB SCRIPTS      
#==================================================
#--check pair-end sequence's id in order or not-------
sub pair_name{
	my $id = shift;
	if ($id =~ /2:N:0:\w+/) {
		 $id=~ s/2:N:0:/1:N:0:/;
		 return $id;
	}else {
		$id=~ s/1:N:0:/2:N:0:/;
		return $id;
	}
}
#---------cut seq from fastq stuff--------------
sub seq {
	my $seq = (split /\n/,shift)[0];
	return $seq;
}

#---------reverse and complement fasta format-----------
sub revcom {
	shift;
	chomp;
	tr/NATCG/NTAGC/;
	$_ = reverse($_);
	return $_;
}
#----------reverse and complement fastq format-------
sub revfq {
	my @aa = split /\n/,shift;
	my $seq = $aa[0];
	my $quality = $aa[2];
	$quality = reverse $quality;
	$seq = &seq ($seq);
	my $new = $seq."\n"."+"."\n"."$quality";
#	print "$new\n";
	return $new;
}

sub output {
	my $str = shift;
	my ($type,$num);
	if ($str =~ /(For|Rev)(\d+)/) {$type = $1 ;$num = $2 ;}
	$type =~ s/rev//;
	my $outfile = "$outdir/"."$num/"."$type"."$num";
}

