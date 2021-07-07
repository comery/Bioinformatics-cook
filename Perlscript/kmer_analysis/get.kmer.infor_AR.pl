#!/usr/bin/perl -w
#use strict;
use File::Basename;
die "Usage;<KMER_log> <KMER_freq> <OUT_dir> <xlim> <ylim> <Title>!" unless @ARGV >5;
my $kmerlog=shift;
my $kmerfreq=shift;
my $outdir=shift;
my $xlim=shift;
my $ylim=shift;
my $title=shift;

#$outdir ||= ""./;
mkdir $outdir unless (-d "$outdir");

open IN,"$kmerlog";# the log from kmer analysis;
open INFOR,">$outdir/infor.txt";
my @A;my $k;
while(<IN>){
	chomp;
	push(@A,$_);
}
for (my$i=0;$i<=$#A;$i++){
	if($A[$i]=~/genome_size/){
		$k=$i;
		last;
	}
}

my @B=split /\s+/,$A[$k];
my @C=split /\s+/,$A[$k+1];
for($i=0;$i<=$#B;$i++){
	print INFOR "$B[$i]\t$C[$i]\n";
}

open IN2,"$kmerfreq";#kmer.freq
my @freq;$i=0;
while(<IN2>){
	chomp;
	my @tmp=split /\s+/;
	$freq[$i]=[@tmp];$i++;
#	push(@freq,@tmp);
}
#for($i=0;$i<=$#freq;$i++){
#	print "$freq[$i][0]\n";
#}
my $kmeruniq=($C[-1]-$freq[0][1]-$freq[1][1]-$freq[2][1]);
print INFOR "kmer_uniq\t$kmeruniq\n";

##########DATA_FOR_PICTURE################################################
open DATA,">$outdir/data.list";
my$freq_acc;my$pro;my$pro_acc;
for($i=0;$i<=$#freq;$i++){
	if($freq[$i][0]==255){
		$pro=$freq[$i][3]/$C[1]*100;
		print DATA "$freq[$i][0]\t";
		print DATA $freq[$i-1][2];
		print DATA "\t100\t$pro\t100\n";
	}else{
		$freq_acc+=$freq[$i][2];
		$pro=$freq[$i][3]/$C[1]*100;
		$pro_acc+=$pro;
		print DATA "$freq[$i][0]\t";
		print DATA $freq[$i][2];
		print DATA "\t".$freq_acc;
		print DATA "\t$pro\t$pro_acc\n";
	}
#	$freq_acc+=$freq[$i][2];
#	$pro=$freq[$i][3]/$C[1]*100;
#	$pro_acc+=$pro;
#	print DATA "$freq[$i][0]\t$freq[$i][2]\t$freq_acc\t$pro\t$pro_acc\n";


}

###############################################################################
open Hr,">$outdir/kmer.draw.R";
print Hr "pdf(file='$outdir/infor.$title.pdf',width=16,height=10)\n";
#print Hr "png(file='$outdir/infor.out.png',width=1200,height=800)\n";
print Hr "layout(matrix(c(1,2,3,4),nr=2))\n";
print Hr "freq <- read.table('$outdir/data.list')\n";
print Hr "par(mai=c(0.7,1.2,0.7,1),font=2,font.lab=2,font.axis=2,lwd=2)\n";
print Hr "xrg <- c(0,$xlim)\n";
#picture1

print Hr "plot(freq[1:$xlim,1],freq[1:$xlim,2],type='l',ylab='Frequency(%)',xlab='Depth',ylim=c(0,$ylim),xlim=c(0,$xlim),lwd=3,col='blue',cex=10,ps=10,cex.main=60)\n";
print Hr "par(new=T,ann=F)\n";
print Hr "plot(freq[1:$xlim,1],freq[1:$xlim,4],type='l',yaxt='n',,lwd=3,col='red',xlim=c(0,$xlim),ylim=c(0,$ylim))\n";
print Hr "axis(4,at=($ylim/2),labels=c('Product(%)'),padj=2.4,tick=FALSE)\n";
print Hr "axis(4, col = 'black', lwd = 1)\n";
print Hr "mtext('$title',line=1.5,side=3)\n";
print Hr "legend('topright',legend=c('Frequency','Product'),col=c('blue','red'),lty=c(1,1))\n";
#picture2
print Hr "par(mai=c(1,1.2,0.7,1),font=2,font.lab=2,font.axis=2,lwd=2)\n";
print Hr "plot(freq[1:$xlim,1],freq[1:$xlim,2],type='l',ylim=c(0,$ylim),xlim=c(0,$xlim),lwd=3,col='blue',cex=0.5)\n";
print Hr "par(new=T,ann=F)\n";
print Hr "plot(freq[1:$xlim,1],freq[1:$xlim,3],type='l',yaxt='n',,lwd=3,col='green',xlim=c(0,$xlim),ylim=c(0,100))\n";
print Hr "axis(4,at=(50),labels=c('Accumulative(%)'),padj=2.4,tick=FALSE)\n";
print Hr "axis(4, col = 'black', lwd = 1)\n";
print Hr "mtext('Depth',side=1,line=3)\n";
print Hr "mtext('Frequency(%)',side=2,line=3)\n";
#table

print Hr "par(mai=c(0.5,0.3,0.5,0.2),font=2,font.lab=2,font.axis=2,lwd=2)\n";
print Hr "plot(1:10,axes=FALSE,pch='',ann=FALSE)\n";
print Hr "abline(h=6)\n";
print Hr "abline(h=6.6)\n";
print Hr "abline(h=5.2)\n";
$k=1;
for($i=0;$i<=6;$i++){
        print Hr "text($k,6.3,'$B[$i]',cex=1.1)\n";
        $k+=1.25;
}
print Hr "text(9.75,6.3,'uniqkmer',cex=1.1)\n";

$k=1;
for($i=0;$i<=6;$i++){
	print Hr "q<-format($C[$i],big.mark=',',scientific=F)\n";
	print Hr "text($k,5.6,q,cex=0.9)\n";
	$k+=1.25;
}
print Hr "q<-format($kmeruniq,big.mark=',',scientific=F)\n";
print Hr "text(9.75,5.6,q,cex=0.9)\n";

#picture3
print Hr "par(mai=c(1,1,0.7,1),font=2,font.lab=2,font.axis=2,lwd=2)\n";
print Hr "plot(freq[1:$xlim,1],freq[1:$xlim,4],type='l',ylab='Product(%)',xlab='Depth',ylim=c(0,$ylim),xlim=c(0,$xlim),lwd=3,col='red',cex=0.5)\n";
print Hr "par(new=T,ann=F)\n";
print Hr "plot(freq[1:$xlim,1],freq[1:$xlim,5],type='l',yaxt='n',,lwd=3,col='green',xlim=c(0,$xlim),ylim=c(0,100))\n";
print Hr "axis(4,at=(50),labels=c('Accumulative(%)'),padj=2.4,tick=FALSE)\n";
print Hr "axis(4, col = 'black', lwd = 1)\n";
print Hr "mtext('Depth',side=1,line=3)\n";
print Hr "mtext('Product(%)',side=2,line=3)\n";
print Hr "\n";

print Hr "dev.off()\n";



system "/share/app/R-3.2.1/bin/R < $outdir/kmer.draw.R --no-save\n";

