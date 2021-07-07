cat common.id1.list|while read name
do
	perl fa_len_cutoff.pl outdir/$name/$name.fa >outdir/$name/$name.del.fa
	muscle -in outdir/$name/$name.del.fa -out outdir/$name/$name.clw -clw -quiet  
	perl ./link_muscle2.pl outdir/$name/$name.clw >outdir/$name/$name.link
	perl change_muscle2fasta.pl outdir/$name/$name.clw >outdir/$name/$name.aln
done
