# get codon position in protein sequence
#python3 phase_PSG.py 13995gene.marmoset.pep psg_region_domain.list > psg.sites
#/hwfssz1/ST_DIVERSITY/PUB/USER/zhouyang/bin/convert/gff2pos.condon.pl
#perl gff2pos.condon.pl 13995gene.marmoset.gff > gff2pos.bed
#python3 pickPos.py gff2pos.bed psg.sites > psg.site.info
#paste psg.site.info scaf.list >tmp
#mv tmp psg.site.info

#
cat region.list|while read a b c
do
	samtools view -h /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/Mapping_10X/1.bwa_mat_curated/merge2one/MatMap10x.merge.sorted.bam $a:${b}-${c} |samtools view -bS - >${a}-${b}-${c}.bam
	perl /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/2.genetic_divergency/SNP_INDEL/bin/Alignmentout_Sum_yct.pl -in ${a}-${b}-${c}.bam -ref /hwfssz5/ST_DIVERSITY/P18Z10200N0160_Hyena/USER/yangchentao/marmoset/0.assembly/curated/mat.genome.curated.fa -type bam -filter no -outpre ${a}-${b}-${c} -outdir ./ -iTools /hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/Bio_packages/iTools_Code/iTools 
done
cat *.xls >all.depth
#rm *.coverage *.pileup *.xls
python3 pickDepth.py all.depth psg.site.info >psg.site.info.depth
