#!/usr/bin/perl -w
use strict;
use Getopt::Long;
=head1 Description

        --sam   "sam file"
        --depthout   "the filename of the depth output"
	--statout "the filename of the summary output"
	--cvgout "the filename of querys' coverage"
        --help  "print out this message"

=cut

my ($Sam,$dOut,$sOut,$cOut,$Help);
GetOptions(
        "sam:s"=>\$Sam,
        "depthout:s"=>\$dOut,
	"statout:s"=>\$sOut,
	"cvgout:s"=>\$cOut,
        "help:s"=>\$Help
);
die `pod2text $0` if ($Help || !defined $cOut || !defined $Sam || !defined $dOut || !defined $sOut);

open IN, "$Sam" || die $!;
my (%hash, %read, %reflen);
while (<IN>){
        chomp;
        if (/^\@SQ/){
                my @a=split;
                $a[1]=~s/SN://;
                $a[2]=~s/LN://;
		$reflen{$a[1]}=$a[2];
		my @s=split /_/, $a[1];
                my @b;
                for (my $i=0; $i<=$a[2]-1; $i++){
                        $b[$i]=0;
                }
                $hash{$a[1]}= [@b];
                $read{$a[1]}=0;
        }
        next if (/^\@/);
        my @c=split;
        next if (($c[1] & 0x0004) or ($c[1] & 0x0200) or ($c[1] & 0x0400));
        my $len=length $c[9];
        my $tmp_j = 1;
        for (my $ii=0; $ii<=$len-1; $ii++){
                my $as=$c[3]-1+$ii;
                if (defined $hash{$c[2]}[$as]){
                        $hash{$c[2]}[$as]++;
                        if ($tmp_j){
						        $read{$c[2]}++;
								$tmp_j = 0;
                        }
                }else{
						print "$c[2]\t$c[0]\n";
                }
        }
}
close IN;

open OUT, ">$dOut" || die $!;
open OUT1, ">$sOut" || die $!;
open OUT2, ">$cOut" || die $!;
my (%cover, %dep);
for my $key (keys %hash){
        print OUT "$key\t$read{$key}".'_reads'."\t";
		print OUT1 "$key\t$read{$key}".'_reads'."\t";
		my $local;
		my $coverlen=0;
		$dep{$key} = 0;
        for my $i (0..$#{$hash{$key}}){
                print OUT "$hash{$key}[$i]\t";
		$dep{$key} += $hash{$key}[$i]; 
				$coverlen++ if ($hash{$key}[$i] >= 1);
				if ($i==0) {
					if ($hash{$key}[0] >= 1) {
						print OUT1 "1";
						if ($hash{$key}[1] >= 1) {
							print OUT1 "-";
						}else{
							print OUT1 ";";
						}
					}
				}else{
					if ($hash{$key}[$i] >= 1) {
						if ($hash{$key}[$i-1] == 0) {
							$local = $i+1;
							print OUT1 "$local";
							if ($hash{$key}[$i+1] >= 1) {
								print OUT1 "-";
							}else{
								print OUT1 ";";
							}
						}elsif ($hash{$key}[$i+1] ==0) {
							$local=$i+1;
							print OUT1 "$local;";
						}
					}
				}
        }
	$dep{$key} = $dep{$key} / $#{$hash{$key}};
	$dep{$key} = sprintf("%.1f", $dep{$key});
        print OUT "\naverage_depth: $dep{$key}\n";
		print OUT1 "\n";
		my @scaf=split /_/, $key;
		if (exists $cover{$scaf[0]}) {
			$cover{$scaf[0]}+=$coverlen;
		}else{
			$cover{$scaf[0]}=$coverlen;
		}
		if (exists $read{$scaf[0]}) {
			$read{$scaf[0]}+=$read{$key};
		}else{
			$read{$scaf[0]}=$read{$key};
		}
}
foreach  my $it (keys %cover) {
		my $coverpercent=$cover{$it}/$reflen{$it}*100;
		my $coverpercent1= sprintf"%0.1f", $coverpercent;
		$coverpercent1 .= '%';
		print OUT2 "$it\t$read{$it}\t$coverpercent1\n";
}
close OUT;
close OUT1;
close OUT2;
print "all done!\n";
