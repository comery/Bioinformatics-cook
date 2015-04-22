#!usr/bin/perl -w 
use strict;
open IN,shift;#the file format is fasta;
open IN1,shift;#the file consist of the id of seq and fam_id queery start ,end;
open IN2,"com_all_best_to_best.id";
my (%aa,%bb,@yy,$fam_id,$id);
$/=">";<IN>;$/="\n";
while (<IN>) {
	chomp;
	$id=$_;
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$/="\n";
	$aa{$id}=$seq;
	#print "$seq\n";
}

while (<IN1>) {
	chomp;
	my @yy=split /\s+/,$_;
	my $key=$yy[1];
	my @kk=splice(@yy,1,1);
	my $xx=join " ",@yy;
	$bb{$key}=$xx;

}

while (my $fam_id=<IN2>) {
	chomp $fam_id;
	my $val=(split /\s/,$bb{$fam_id})[0];
	my $start=(split /\s/,$bb{$fam_id})[1];
	my $end=(split /\s/,$bb{$fam_id})[2];
	my $fa=$aa{$val};
	#print "$fa";
	if($end > $start){
		my $len=$end-$start+1 ;
    	my $querry=substr($fa,$start-1,$len);
    	print  ">$val\_$fam_id\n$querry\n";
		}else{
		my $len=$start-$end+1;
		my $querry=substr($fa,$end-1,$len);
		print ">$val\_$fam_id\n$querry\n";
		}

}
