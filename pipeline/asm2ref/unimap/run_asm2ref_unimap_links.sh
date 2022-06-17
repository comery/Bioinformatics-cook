if [[ $# != 3 ]] ; then
	echo "Usage : $0 ref.fasta query.fasta outprex"
	exit 1
fi

ref=$1
qry=$2
outpre=$3
src=$(cd $(dirname $0);pwd)
# unimap
/hwfsxx1/ST_HN/P18Z10200N0197/yangchentao/software/T2T/unimap/unimap -cxasm5 --cs -t8 $ref $qry > $outpre.unimap.paf
# generate link.txt and karyotype.txt
python3 $src/preCircosLinkFromPaf.sort.py -paf $outpre.unimap.paf -filter 
# get conf 
cp $src/confs/run_asm2ref_unimap_links.conf .
# plot cirocs
circos -conf run_asm2ref_unimap_links.conf
