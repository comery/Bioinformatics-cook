![example](http://ogj9k5cjf.bkt.clouddn.com/NC_031379.png)
#### summary
circos is very powerful software to illustrate genetic information, especially for mitochondrion and chloroplast genomes. Here I provide a useful perl script to convert genome annotation with genebank format to all features as circos input, and then draw a circos image.

#### requirement

- overall, circos must be installed in your computer.

actually, there are two versions, one is “standard”, and the other is “auto_depth”.
- for “standard” mode, you just need a genebank input file, and general configure file, like below.
```
# where circos you install
circos_path		=	/usr/local/bin/circos
bwa             =   /usr/local/bin/bwa
samtools        =   /usr/local/bin/samtools
#color
cds				=   102,194,165
rRNA			=	252,141,98
tRNA			=	123,50,148
# whether draw GC content circle
gc				=	yes
win				=	50
gc_fill			=	146,197,222
# whether draw depth abundance circle
depth			=	yes
depth_fill		= 	5,113,176
depth_file		=	../example.depth.txt
# whether draw base around circle
base 			=	no
# locus name's clolor showed on center of circle
locus_color		=	black
# gene name label color
label_color		=	black
# image
outdir			=	./example
png				=	yes
svg				=	yes
# color or file
background		=	white
#background		=	./background.png
```

- for “auto_depth” mode, you also need a genebank input file and configure file, but what you should notice carefully is path setted in head of configure file - "circos_path", "bwa", and "samtools".

```
# where circos you install
circos_path		=	/usr/local/bin/circos
bwa             =   /usr/local/bin/bwa
samtools        =   /usr/local/bin/samtools
#color
cds				=	102,194,165
rRNA			=	252,141,98
tRNA			=	123,50,148
# whether draw GC content circle
gc				=	yes
win				=	50
gc_fill			=	146,197,222
# whether draw depth abundance circle
depth			=	no
depth_fill		= 	5,113,176
fq				=	test1.fq,test2.fq

# whether draw base around circle
base 			=	no
# locus name's clolor showed on center of circle
locus_color		=	black
# gene name label color
label_color		=	black
# image
outdir			=	./outdir
png				=	yes
svg				=	yes
# color or file
background		=	white
#background		=	./background.png
```
as a result, paired fastq file should give in before file.


now, you can enjoy it.

```
perl draw_circos_for_mitogenome_auto_depth.pl -gb ../example.gb -conf mitogenome.auto_depth.conf 

perl draw_circos_for_mitogenome_standard.pl -gb ../example.gb -conf mitogenome.standard.conf

```

