#!usr/bin/perl -w
use strict;

die "usage: perl $0 <.gb> <output filename>\n" unless(@ARGV == 2);

open (GB, $ARGV[0]) or die "cannot open $ARGV[0] $!\n";
open (TFA, ">$ARGV[1].fa");
open (PTFA, ">$ARGV[1].PT.fa");
open (AAS, ">$ARGV[1].aas");
open (RRNA, ">$ARGV[1].rRNA.fa");
open (TRNA, ">$ARGV[1].tRNA.fa");
open AA1, ">$ARGV[1].gb.erro1";
open AA2, ">$ARGV[1].sub.split.erro";


$/ = "//\n";

my ($nc_id,$sp,$seq,$str,@Unit,$direction,$unit,$start,$start2,$end,$len,@split_gene,$trna,$sub_seq,$rrna,$cds,$gi,$aaa_seq,$aaa_len);

while (<GB>) {
	($nc_id) = $_ =~ /ACCESSION\s+(\w+)\n/;
	($sp) = $_ =~ /ORGANISM\s+(.*)\n/;
	$sp =~ s/\s+/_/g;
	($seq) = $_ =~ /ORIGIN\s*\n([a-z0-9\s]+)\s*/;
	$seq =~ s/\s*//g; 
	$seq =~ s/\d+//g;
	print TFA ">$nc_id\_mitogenome|$nc_id|$sp-mitogenome\n$seq\n";
	($str) = $_ =~ /(source[a-zA-Z0-9"-_\(\)\s=\.:\/]+)ORIGIN\s*\n[a-z0-9\s]+/;
	unless ($str){print AA1 "$str*\n\n*"};
	@Unit = split /\s*tRNA\s* | \s*rRNA\s* | \s+?CDS\s+?/, $str;
	shift @Unit;
    foreach(@Unit) {
        $unit = $_;
		#print "UUUUUUUUUUU\n$unit\n UUUUUUUUUUU\n";
		if (/complement\(/){
			$direction = 0;
		}
		else{
			$direction = 1;
		}
        if (/tRNA/){
			($start) = $unit =~ m/(\d+)\.\./;
			$start = $start - 1;
			$start2 = $start + 1;
			($end) = $unit =~ m/\.\.(\d+)/;
			$len = $end - $start;
			@split_gene = split /\s+gene\s+/, $unit;
			($trna) = $split_gene[0] =~ m/\/product="(.*)"\s*/;
			$sub_seq = substr ($seq,$start,$len);
			print TRNA ">$nc_id\_$trna|$nc_id|$trna|$sp|0|$nc_id\_mitogenome|$start2|$end|$direction\n$sub_seq\n";
			$start = undef;$end = undef;$len = undef;@split_gene=();$trna=undef;$sub_seq = undef;
		}
		elsif (/translation="/){
			($start) = $unit =~ m/(\d+)\.\./;
			$start = $start - 1;
			$start2 = $start + 1;
			($end) = $unit =~ m/\.\.(\d+)/;
			$len = $end - $start;
			@split_gene = split /\s+gene\s+/, $unit;
			($cds) = $split_gene[0] =~ m/\/product="(.*)"\s*/;
			$cds = &CDS($cds);
			$cds =~ s/NDL/ND4L/g;
			($gi) = $split_gene[0] =~ m/\/db_xref="GI:(\d+)"\s*/;
			$sub_seq = substr ($seq,$start,$len);
			($aaa_seq) = $split_gene[0] =~ m/\/translation="([A-Z\s]+)"/;
			$aaa_seq =~ s/\s+//g;
			$aaa_len = length ($aaa_seq);
			print PTFA ">gi_$gi\_$cds|$nc_id|$cds|$sp|0|$nc_id\_mitogenome|$start2|$end|$direction\n$sub_seq\n";
			print AAS ">gi_$gi\_$cds\_$sp\_$aaa_len\_aa\n$aaa_seq\n";
			$start = undef;$end = undef;$len = undef;@split_gene=();$cds=undef;$gi=undef;$sub_seq=undef;$aaa_seq=undef;$aaa_len=undef;
		}
		elsif (/r.*?RNA/){
			($start) = $unit =~ m/(\d+)\.\./;
			$start = $start - 1;
			$start2 = $start + 1;
			($end) = $unit =~ m/\.\.(\d+)/;
			$len = $end - $start;
			@split_gene = split /\s+gene\s+/, $unit;
			($rrna) = $split_gene[0] =~ m/\/product="(.*)"\s*/;
			$rrna = &rRNA($rrna);
			$sub_seq = substr ($seq,$start,$len);
			print RRNA ">$nc_id\_$rrna|$nc_id|$rrna|$sp|0|$nc_id\_mitogenome|$start2|$end|$direction\n$sub_seq\n";
			$start = undef;$end = undef;$len = undef;@split_gene=();$rrna=undef;$sub_seq = undef;
		}
		elsif (/\d+S/){
			($start) = $unit =~ m/(\d+)\.\./;
			$start = $start - 1;
			$start2 = $start + 1;
			($end) = $unit =~ m/\.\.(\d+)/;
			$len = $end - $start;
			@split_gene = split /\s+gene\s+/, $unit;
			($rrna) = $split_gene[0] =~ m/\/product="(.*)"\s*/;
			$rrna = &rRNA($rrna);
			$sub_seq = substr ($seq,$start,$len);
			print RRNA ">$nc_id\_$rrna|$nc_id|$rrna|$sp|0|$nc_id\_mitogenome|$start2|$end|$direction\n$sub_seq\n";
			$start = undef;$end = undef;$len = undef;@split_gene=();$rrna=undef;$sub_seq = undef;
		}
		else {
			print AA2 "*****erro*****\n$unit******erro*****\n";
		}
	}
	$nc_id=undef;$seq=undef;$str=undef;@Unit=();
}

close TFA;
close GB;
close PTFA;
close AAS;
close TRNA;
close RRNA;
close AA1;
close AA2;

my $gene;
sub CDS {
	$cds = $_[0];
	
	if ($cds =~ /ATP.*?\s(8)/ || $cds =~ /ATP.*?\s(6)/){
		$cds = join "",'ATP',$1;
	}
	elsif ($cds =~ /NADH dehydrogenase .*([123456L]+)/){
		$cds = join "",'ND',$1;
	}
	elsif ($cds =~ /adenosine .* (\d)/){
		$cds = join "",'ATP',$1;
	}
	elsif ($cds eq 'cytochrome oxidase subunit I' || $cds eq 'cytochrome c oxydase subunit I' || $cds eq 'cytochrome c oxidase subunit I' || $cds =~ /cytochrome.* 1/){
		$cds = 'COX1';
	}
	elsif ($cds eq 'cytochrome oxidase subunit II' || $cds eq 'cytochrome c oxydase subunit II' || $cds eq 'cytochrome c oxidase subunit II' || $cds =~ /cytochrome.* 2/){
		$cds = 'COX2';
	}
	elsif ($cds eq 'cytochrome oxidase subunit III' || $cds eq 'cytochrome c oxydase subunit III' || $cds eq 'cytochrome c oxidase subunit III' || $cds =~ /cytochrome.* 3/){
		$cds = 'COX3';
	}
	elsif ($cds =~ /cytochrome b/){
		$cds = 'CYTB';
	}
	$cds;
}

sub rRNA {
	$gene = $_[0];
	if ($gene =~ /s-rRNA/ || $gene =~ /12S/i){
		$gene = '12SrRNA';
	}
	elsif ($gene =~ /l-rRNA/ || $gene =~ /16S/i){
		$gene = '16SrRNA';
	}
	$gene;
}
