# fastaKit

```text
usage: fastaKit [-h] [-n <STR>] [-ex <STR>] [-rg <STR>] [-nl FILE] [-exl FILE] [-rgl FILE] [-ap <STR>] [-pos <INT>] [-lgt <INT>] [-lle <INT>] [-a2u] [-rw]
                [-tr] [-codon CODON] [-cutf <INT>] [-cuts <INT>] [-o <STR>] [-z] [-sta]
                <File>

Description

    Fasta Kit for dealing with fasta file, including name modification,
    filtering sequences by length, format changes, subfiles retrevied, etc.
    I include all kinds of tools as possible as I can.

Usage
    # select sequences which contain 'HOMO' in ID
    python3 fastaKit  -n HOMO -o output.fa test.fa
    # add a name 'prefix' to name ID for each sequence
    python3 fastaKit  -ap HOMO -o output.fa test.fa
    # filtering sequence by length (>= 100bp)
    python3 fastaKit  -lgt 100 -o output.fa test.fa
    # filtering sequence by length (>= 100 and <= 150)
    python3 fastaKit  -lgt 100 -lle 150 -o output.fa test.fa
    # aligned sequence to unaligned
    python3 fastaKit  -a2u  -o output.fa test.fa
    # cut sequence into [INT] files
    python3 fastaKit  -cutf 10 -o outdir test.fa

positional arguments:
  <File>        input fasta file

optional arguments:
  -h, --help    show this help message and exit
  -n <STR>      target ID exactally
  -ex <STR>     excluded name
  -rg <STR>     target regx word to find sequences
  -nl FILE      selecting sequences by a ID list
  -exl FILE     selecting sequences by an excluding list
  -rgl FILE     selecting sequences by a regular expression list
  -ap <STR>     add a specific prefix to ID
  -pos <INT>    position of tag you added
  -lgt <INT>    deal with by length, select sequences by larger than [INT] bp
  -lle <INT>    deal with by length, select sequences by smaller than [INT] bp
  -a2u          change format, aligned sequence to unaligned
  -rw           change format, multi-line to single line
  -tr           translate fasta into protein (CDS only)
  -codon CODON  codon table usage (CDS only, with -tr)
  -cutf <INT>   cut fasta into files by file number
  -cuts <INT>   cut fasta into files by specific sequence number
  -o <STR>      output file or outdir for -cutf/-cuts
  -z            output a gzip type file
  -sta          basic statistic of fasta
```