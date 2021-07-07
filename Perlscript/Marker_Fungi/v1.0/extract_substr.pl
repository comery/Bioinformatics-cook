 ######################################
 #This perl is to extract the seq from
 # the blast.out ,For your need ,you can 
 #change the number in "$a[17]"and "$a[18]"
 ######################################
 #!  usr/bin/perl -w
 use strict;
 open IN,shift;
 open IN1,shift;
 my (@a,%hash,$start,$end);
 $/=">";<IN>;
 while (<IN>){
 	chomp;
 	@a=split /\n/,$_;
 	$hash{$a[0]}=$a[1];
 }
local $/="\n";
while (<IN1>) {
	@b=split(/\s+/,$_);
	my $id=$b[0];
	my $seq=$hash{$id};
    my $len=$end>$start ?$end-$start+1:$start-$end+1;
    my $querry=substr($seq,$start-1,$len);
    print  ">$id\n$querry\n";
  
}
close IN;
close IN1;
