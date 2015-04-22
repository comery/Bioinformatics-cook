#!/usr/bin/perl
use strict;
use warnings;
die "Usage: perl 1(gene family num)" if(@ARGV!=1);
my $all_cds_file="/ifs4/NGB_ENV/USER/yangchentao/Fungi-genomes/Verticillium/Verticillium_case/ortholog/ortholog.rename.fa";
my $single_copy_file="/ifs4/NGB_ENV/USER/yangchentao/Fungi-genomes/Verticillium/Verticillium_case/ortholog/id.sort.list";
my $gene_family=shift;
`mkdir -p ./cluster/$gene_family`; 

`less $single_copy_file | awk '\$1=="$gene_family"' | cut -f 2-125 | sed 's/ /\\n/g' > ./cluster/$gene_family/id_list`;

`perl /ifs4/NGB_ENV/USER/yangchentao/software/fishInWinter.pl  -bf table -ff fasta ./cluster/$gene_family/id_list $all_cds_file > ./cluster/$gene_family/id.cds.fa`;

`/ifs4/NGB_ENV/USER/yangchentao/software/muscle -in  ./cluster/$gene_family/id.cds.fa  -msf  -out ./cluster/$gene_family/id.cds.fa.out`;

` perl /ifs4/NGB_ENV/USER/yangchentao/Fungi-genomes/Verticillium/Verticillium_case/data/hmm/change.pl ./cluster/$gene_family/id.cds.fa.out > ./cluster/$gene_family/id.cds.fa.out.sto`;

`/ifs5/PC_PA_UN/ENV/USER/zhoulili/software/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild ./cluster/$gene_family/id.cds.fa.out.hmm ./cluster/$gene_family/id.cds.fa.out.sto`;


