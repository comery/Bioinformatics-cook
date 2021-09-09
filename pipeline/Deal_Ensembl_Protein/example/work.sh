date
bin="../bin"
genome="ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf="ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz"
pep="ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"

# download above three files at first
pre="homo_sapiens"
perl $bin/renameGenome.pl $genome >$pre.fa
perl $bin/gtf2gff.pl $gtf >$pre.all.gff
perl $bin/check_orf_for_gff.pl $pre.all.gff $pre.fa >$pre.all.gff.orf
perl $bin/bestProteinFromEnsembl_87.pl $pep $pre.all.gff.orf 1 >$pre.pep.best
perl $bin/select_gff.pl $pre.pep.best $pre.all.gff >$pre.best.gff
perl $bin/getGene.pl $pre.best.gff $pre.fa >$pre.best.cds
perl $bin/check_orf_for_cds.pl $pre.best.cds | awk '$4>0' >$pre.best.cds.preStop
perl $bin/cds2aa.pl $pre.best.cds >$pre.best.pep
perl $bin/select_orf.pl $pre.best.pep $pre.all.gff.orf >$pre.best.gff.orf
awk '$8==1 && $9==1' $pre.best.gff.orf >$pre.best.gff.orf.intact
rm $pre.all.gff $pre.all.gff.orf $pre.pep.best
date
