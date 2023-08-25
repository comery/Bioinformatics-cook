this script is designed for downstream when generated the results from PAML 2.0 pipeline.
input.list must contain the list of family dirtory, for example,
$ cat input.list
family00306

family00306 is directory name, and it contains two necessary file:
- family00306.seq
- H1.mlc

it will draw the phylogenetic tree and specific region in alignment which you defined by parameters.
Meanwhile, it will output family00306.tree, family00306.fasta two files.

