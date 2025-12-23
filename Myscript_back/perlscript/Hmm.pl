#!/usr/bin/perl

use strict;
use warnings;

die "Usage: perl 1(gene family num)" if(@ARGV!=1);
my $all_cds_file="../orthologs.rname.fas";
my $single_copy_file="../cluster_65_80_1716.nr";
my $gene_family=shift;
`mkdir ../$gene_family`; 
`less $single_copy_file | awk '\$1=="$gene_family"' | cut -f 2-125 | sed 's/ /\\n/g' > ../$gene_family/id_list`;
`perl /ifs5/PC_PA_UN/ENV/USER/zhoulili/bin/command/fishInWinter.pl  -bf table -ff fasta ../$gene_family/id_list $all_cds_file > ../$gene_family/id.cds.fa`;
`/ifs1/ST_SINGLECELL/USER/jiangrunze/bin/muscle -in  ../$gene_family/id.cds.fa  -msf  -out ../$gene_family/id.cds.fa.out`;
` perl /ifs5/PC_PA_UN/ENV/USER/zhoulili/1KITE/venom/Hym-HMM/Hym/HMM/change.pl ../$gene_family/id.cds.fa.out > ../$gene_family/id.cds.fa.out.sto`;
`/ifs5/PC_PA_UN/ENV/USER/zhoulili/software/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild ../$gene_family/id.cds.fa.out.hmm ../$gene_family/id.cds.fa.out.sto`;
` perl /ifs5/PC_PA_UN/ENV/USER/zhoulili/software/hmmer-3.1b1-linux-intel-x86_64/translate.pl -infile /ifs5/PC_PA_UN/ENV/USER/zhoulili/1KITE/AQIS-Fungi/Fungi-genomes/test/Verticillium_dahliae.rname.fas  -trunc 0  -outfile Verticillium_dahliae.out`;
`/ifs5/PC_PA_UN/ENV/USER/zhoulili/software/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmsearch --domtblout ../$gene_family/Verticillium_dahliae.result   ../$gene_family/id.cds.fa.out.hmm  /ifs5/PC_PA_UN/ENV/USER/zhoulili/1KITE/AQIS-Fungi/Fungi-genomes/test/Verticillium_dahliae.out`;
