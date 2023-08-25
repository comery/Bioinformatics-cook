#!/usr/bin/perl

=head1 cdsgff2seq.pl

read GFF, genome fasta; write CDS sequence, aa translation (with check for exon phase).
write gene/exon offset subranges.

=head1 WARNING Bio::DB::Fasta needs even lines

** hidden problem in switch to Bio::DB::Fasta for indexing:
it requires fasta have all same line width, and older (>2007) versions
dont report when this isn't true, but return wrong sequence.

Reformat to even width fasta if needed:

  cat old.fasta | perl -ne\
  'if(/^>/){ if($a){ $a =~ s/(.{1,50})/$1\n/g; print $a,"\n"; } $a=""; print; }\
  else{ chomp; $a.=$_;} END{if($a){$a =~ s/(.{1,50})/$1\n/g; print $a,"\n"; }}'\
  > new.fasta


=head1 usage

    -- CDS sequence
  cdsgff2seq -a dna -t exon region.gff chromfa.dir/ > out.cdna
    -- Amino translation
  cdsgff2seq -a aa -t CDS region.gff chromfa.dir/ > out.pep

    -- introns of CDS sequence
  cdsgff2seq.pl -a intron -t CDS -gff dpulex-sc4-gno.gff dpulex-sc4.fa > dpulex-sc4-intron.fa


    -- Subsequence ranges around gene,exon boundaries
  cdsgff2seq.pl -a genestart -offset=-99,9 -t exon -gff ${dp}_scaf7*.gff -fasta $sc/${dp}/perchr/ > ${dp}-sc7g5.fa
  cdsgff2seq.pl -a intronstart -offset=-19,20 -t exon -gff ${dp}_scaf7*.gff -fasta $sc/${dp}/perchr/ > ${dp}-sc7in1.fa
  cdsgff2seq.pl -a intronend -offset=-19,20 -t exon -gff ${dp}_scaf7*.gff -fasta $sc/${dp}/perchr/ >  ${dp}-sc7in2.fa

 options: -v verbose ; -t CDS restrict to this feature (CDS,exon default)
        -a dna or -a aa : option to output aa translations, cds dna

  -- genbank output : fixme to make simpler options?

  cdsgff2genbank.pl -t='CDS,exon,gene' -a=genbank -gff=dgri_s15252c.gff -fasta=dgri_caf060210.fa > dgri_s15252_genes.gb

  
=head1 Requirements
  
  Uses BioPerl Bio::SeqIO to read sequence files for exon subranges
  Inputs: 
  1. GFF v3 for feature locations (exons, CDS).  Parent= or ID= must
     distinguish genes.
  2. Directory of FastA genome sequence matching GFF Reference field
    (e.g. faSplit byname genome.fasta)
     Uses dna fasta file per scaffold/chromosome as $dnapath/$ref.fa
    
    
=head1 Author
  
  d.gilbert, aug2006, gilbertd@indiana.edu
  updated march 2007 for gene/exon subsets
  from cdsphase.pl, august 2006

=head1 subset seqs

 07mar addition: subset sequence output options
 choose range around exons to output: 5prime, 3prime,intron regions
   +/- n bases from boundary
   -a genestart geneend intronstart intronend : fixme
   -offset -50,10
   
  * want option to do several/all -astart/end -offset x,y to separate files 
    in one pass
    
=head1 example subsets

  # get some gff gene models (est6 six-pack)
  curl 'http://insects.eugenes.org/species/cgi-bin/gbrowse/dvir/?name=scaffold_12855:9338335..9363334;label=NCBI_GNO;plugin=GFFDumper;plugin_do=Go' > dvir-est6gno2.gff

  set dp=dvir
  perl cdsgff2seq.pl -a genestart -o=-99,9 -t exon -gff ${dp}-est6gno.gff -fasta $sc/${dp}3/perchr/  > ${dp}-gb1.fa
  perl cdsgff2seq.pl -a geneend -o=-9,99 -t exon -gff ${dp}-est6gno.gff -fasta $sc/${dp}3/perchr/ > ${dp}-ge1.fa

  perl cdsgff2seq.pl -a intronstart -o=-19,20 -t exon -gff ${dp}-est6gno.gff -fasta $sc/${dp}3/perchr/ > ${dp}-in1.fa
  perl cdsgff2seq.pl -a intronend -o=-19,20 -t exon -gff ${dp}-est6gno.gff -fasta $sc/${dp}3/perchr/ > ${dp}-in2.fa

  set dp=dpulex
  perl cdsgff2seq.pl -a genestart -o=-99,9 -t exon -gff ${dp}_scaf7*.gff -fasta $sc/${dp}4/perchr/ > ! ${dp}-sc7g5.fa
  perl cdsgff2seq.pl -a intronend -o=-19,20 -t exon -gff ${dp}_scaf7*.gff -fasta $sc/${dp}4/perchr/ > ! ${dp}-sc7in2.fa
  perl cdsgff2seq.pl -a intronstart -o=-19,20 -t exon -gff ${dp}_scaf7*.gff -fasta $sc/${dp}4/perchr/ > ! ${dp}-sc7in1.fa

  # use WebLogo for intron/exon junction and 5' logos
  http://bespoke.lbl.gov/weblogo/
  
  # or find motifs with ..........
  See http://www.sanger.ac.uk/Software/analysis/nmica/
    for analyzing gene subset sequences for motifs.   
  set mp=geall ; 
  $nm/makemosaicbg -seqs $mp.fa -mosaicClasses 1 -mosaicOrder 1 -out $mp.sbg
  $nm/motiffinder -out $mp-mot.xms -seqs $mp.fa -backgroundModel $mp.sbg -numMotifs 3 

   See also
   Down TA, Bergman CM, Su J, Hubbard TJP (2007) 
   Large-scale discovery of promoter motifs in Drosophila melanogaster. 
   PLoS Comput Biol 3(1): e7. 
   doi:10.1371/journal.pcbi.003000

=cut

use strict;
use lib("/bio/argos/common/perl/lib/");
##use Bio::SeqIO; # replace this with that
use Bio::DB::Fasta;

my $verbose= 0;
my $outh = *STDOUT;
my $gffin= *STDIN;

my $MIN_CDS_LEN = 60;
my $MIN_AA_LEN = 20;
my $trustphase= 0;
my $GENE_EXPAND = 600;

# each predictor uses diff feature types; do we also need mRNA,gene parent?
my @cdsTypes= qw(CDS exon); # FIXME
my %cdsTypes= map { $_ => 1; } @cdsTypes;
my %seenType=();
my $seqtype= "dna"; # dna/cds aa/protein
my $fasta_db= undef;

# expect args:  cdsgff2seq.pl  data.gff[.gz] dnafastapath/
my %args=();

use Getopt::Long;
my $ok= &GetOptions( \%args,
  "gff=s" ,
  "fasta|dnafastapath=s" ,
  "t|types|features=s" , # gff features to keep: CDS, exon, ...
  "a|seqtype=s" , ## dna,aa,protein,...  
                  ## a added: genestart,geneend,intronstart,intronend
                  ## add genbank/genespan ?? syntax
  "o|offset=s" , ## -50..0/genestart ; 0..30/geneend ; -10..40/intronstart ; -40..10/intronend
  "MIN_CDS_LEN=i" , ## -50..0/genestart ; 0..30/geneend ; -10..40/intronstart ; -40..10/intronend
  "casechange!" , ## for -n .. +n, change seq case to reflect -/+ boundary
  "v|verbose!" ,
  "gnomonflags!",
  "pasaflags!",
  "prettyfa!",
  "runfast!",
  "debug!" ,
  );

my $gff= $args{gff} || shift(@ARGV);
my $genomefasta= $args{fasta} || shift(@ARGV);

die "Usage:  cdsgff2seq.pl [-verbose] [-a=genbank|dna|protein|intron] [-t=CDS] -gff=annot.gff[.gz]  -fasta=genome-fasta-file-or-folder/\n
     'perldoc cdsgff2seq.pl' for details\n" 
  unless($ok && $gff && -e $genomefasta);

if($gff =~ /\.gz$/) { $ok= open(GFF, "gunzip -c $gff|"); $gffin = *GFF; }
elsif($gff =~ /stdin/) { $ok=1; $gffin = *STDIN; }
else { $ok= open(GFF,$gff); $gffin = *GFF; }

my $debug= $args{debug};
my $casechange= $args{casechange};
my @offset= ($args{o}) ? split(/[,;.]+/,$args{o}) : ();
my $usegnomonflags= $args{gnomonflags} || 0;
my $usepasaflags=$args{pasaflags} || 0;
my $prettyfa=$args{prettyfa} || 0;
my $runfast= defined $args{runfast} ? $args{runfast} : 1;
$MIN_CDS_LEN= $args{MIN_CDS_LEN} || $MIN_CDS_LEN;

if($args{v}){ $verbose= $args{v}; }
if($args{t}){ @cdsTypes=split(/[,;.-]/,$args{t}); }
if($args{a}){ $seqtype=$args{a}; }
$verbose= 1 if($debug);

warn "# cdsgff2seq: gff=$gff; seq=$genomefasta\n" if $verbose;

processGff($gffin,$genomefasta,$outh);
close($gffin);
close($outh);

#--------------------


sub get_dna {
  my($fasta, $ref, $start, $stop)= @_; #, $fasta_db_ref
  
  unless( $ref && $stop>0) {
    warn "need ref-ID, start, stop location\n"; return;
    }
    
  my $havedb= ref $fasta_db;
  if($havedb and !$runfast) {
    $havedb= $fasta_db->index_name() eq $fasta_db->index_name($fasta,-d $fasta);
    }
  unless($havedb) {
    #? check patched Bio::DB::Fasta for even line width test?
    my $db = eval { Bio::DB::Fasta->new($fasta); }
      or die "$@\nCan't open sequence file(s). "; # and return;
    $fasta_db= $db;  
    $havedb= 1;
    }
  
  my $seq = $fasta_db->seq($ref, $start => $stop) 
      or return; ## die "cant locate seq of $ref:$start-$stop\n";#? and return;
  $seq= $seq->seq if(ref $seq); # is this weird bioperl change here or not
  return $seq;
}


sub processGff {
  my($gffHandle,$genomefasta,$outh)= @_;

  %seenType=();
  %cdsTypes= map { $_ => 1; } @cdsTypes;
  
  my($faIO, $refdna)= (undef, undef);
  my($l_id, $l_ref, $skipref)=("")x3;
  my( @feats, @cds, %idbad, $gffver);
  $gffver=0;
  $gffver=3; # assume ok when? FIXME
  
  # any patch for $ref to gff-ref ?  e.g. dpse 'Ch'; dper/dsec super/scaffold; ...
  ## assume cds grouped by ID in gff

  while(my $gffin= <$gffHandle>) {
    unless($gffin =~ m/^\w/) { 
      if($gffin =~ /##gff-version\s+(\d+)/) { $gffver=$1; }
      next; }  
    
    unless($gffver >= 3) {
      die "This program requires GFF version >=3 with ##gff-version header\n" 
       ."Found gff-version: $gffver\n at $gffin\n" ;
      }
      
    ## need to collect all CDS_exons/gene-mRNA first ...
    my @gff= split "\t",$gffin,9;
    my $ref = $gff[0];
    my $type= $gff[2];
    chomp($gff[8]);
    my $id= $gff[8];  $id =~ s/^(ID|Parent)=//; $id =~ s/[;,\s].*$//;
    $seenType{$type}++;
    # check here for exon AND CDS : drop exon if have CDS
    #?? %cdsTypes=("CDS" => 1) if( $type eq "CDS" ); #$cdsTypes{"exon"}= 0 

    # dang it only on mRNA need for CDS: id-bad flag
    if($usegnomonflags && $gff[8] =~ /flags=([^;]+)/){ 
      my $fl=$1; my $isfull= ($fl=~/Start/i && $fl=~/Stop/i); 
      $idbad{$id}++ unless($isfull);
      }
    if($usepasaflags && $gff[8] =~ /status=([^;]+)/){ 
      my $fl=$1; my $isfull= ($fl=~/complete/i); 
      $idbad{$id}++ unless($isfull);
      }
#     if($dupflags && $gff[8] =~ /flags=([^;]+)/){ 
#       my $fl=$1; my $isdup= ($fl=~/dupl/i); 
#       $idbad{$id}++ if($isdup);
#       }
    next if($idbad{$id});
    
    if($id ne $l_id && @cds) {
      my $cdsnew= cds2seq( $genomefasta, $l_ref, \@cds, $seqtype); # sorted @cds / gene; not all @feats
      # my $cdsnew= cds2seq( $refdna, \@cds, $seqtype); # sorted @cds / gene; not all @feats
      # foreach my $gff (@$cdsnew) { print $outh join("\t",@$gff),"\n"; }
      print $outh $cdsnew if $cdsnew; # fully formatted?
      @cds=(); # @feats=(); 
      }

    if($ref ne $l_ref) {
      $skipref=0;
      ## replace this with Bio::DB::Fasta use ... 
      
#       undef $refdna; undef $faIO;
#       my $fafile= "$genomefasta/$ref.fa";
#       # check here for $ref.fa variations: dpse: Ch2..ChX..
#       unless(-e $fafile) {
#         my $fref= $ref;
#         if($fref =~ /^Ch/) { $fref =~ s/^Ch//; }
#         elsif($fref =~ /^super/) { $fref =~ s/^super/scaffold/; }
#         $fafile= "$genomefasta/$fref.fa"; 
#       }
#       unless(-e $fafile) {
#         warn "*** Cant find dna.fa for $ref at $fafile\n";
#         $skipref=1;  
#       } else {
#         $faIO = Bio::SeqIO->new(-file => $fafile,-format => 'Fasta');
#         $refdna = $faIO->next_seq(); # wait till have gff ?
#       }
      
      
      # warn "# cdsgff2seq[$ref]: $fafile\n" if $verbose;
      $l_ref= $ref;
      }
      
    ## fixme for "intron" type: exon,CDS diff
    
    if( $cdsTypes{$type} && !$skipref) { push(@cds,\@gff); } 
    # else { print $outh $gffin; }   #? dont need to keep
    
    $l_id= $id; $l_ref= $ref;
    }
    
  if(@cds) {
    my $cdsnew= cds2seq( $genomefasta, $l_ref, \@cds, $seqtype); # sorted @cds / gene; not all @feats
    ##my $cdsnew= cds2seq( $refdna, \@cds, $seqtype); # sorted @cds / gene; not all @feats
    # foreach my $gff (@$cdsnew) { print $outh join("\t",@$gff),"\n"; }
    print $outh $cdsnew if $cdsnew; # fully formatted?
    @cds=();  
    }
    
  undef $refdna; undef $faIO;
}


=item togenbank

  Output in genbank format of * gene-centric * features + sequence.
  Now only useful for genespan with mRNA (and all exons) plus CDS (CDS_exons only)
  add possible other gene-centric features from gff ?

  min format:

LOCUS       geneid        38001 bp    DNA             INV       23-Feb-2000
DEFINITION  ... geneid, gff Note, other?
COMMENT     Converted from gff+fasta with cdsgff2genbank version n.

FEATURES             Location/Qualifiers                    
     source          1..38001
                     /organism="Drosophila melanogaster"

     mRNA            join(3635..4200,15849..16693,18104..18242,
                     18372..18917,19077..19284,24304..24446,25225..25374,
                     29379..29542,29661..29824,32733..34093)
                     /comment="... note ... "
                     /label="Ace-RA|mRNA"
                     /symbol="Ace-RA"

     CDS             join(16371..16693,18104..18242,18372..18917,
                     19077..19284,24304..24446,25225..25374,29379..29542,
                     29661..29824,32733..32845)
                     /aa_size=649
                     /derived_from="Ace-RA"
                     /evidence=predicted
                     /label="Ace-P1|CDS"
                     /symbol="Ace-P1"
                     /translation="MAISCRQSRVLPMSLPLPLTIPLPLVLVLSLHLSGVCGVIDRLVV
                     QTSSGPVRGRSVTVQGREVHVYTGIPYAKPPVEDLRFRKPVPAEPWHGVLDATRLSATC
                     LIYICAALRTKRVF"

ORIGIN
        1 acagagcagt cattacaaga atcaaaatgt ctaatctcaa ctttcaagtt ctagtttatt
    38001 a
//

=cut

sub exons2loc {
  my($exonsGFF,$dnastart)= @_;
  my $loc=""; my $iscomp=0;
  foreach my $ex (@$exonsGFF) {
    my($ref,$src,$type,$start,$stop,$score,$strand,$phase,$attr)= @{$ex};
    $loc.="," if($loc);
    ($start,$stop)= ($start-$dnastart, $stop-$dnastart);
    $loc .= "$start..$stop";
    $iscomp=1 if($strand eq "-" or $strand<0);
  }
  if($loc=~/,/) { $loc="join($loc)"; }
  if($iscomp) { $loc= "complement($loc)"; } #?? does aug want both?
  return $loc;
}

sub wrapline {
  my( $v, $tab)= @_;
  if((my $len= length($v)) > 50) { # 50  # still problems w/ augustus reading these locs
    my $i= 0; my $w="";
    while($i <  $len) {
      my $b= 49;   
      if( $b + $i >= $len) { 
        $w .= substr($v,$i,$len - $i);
        $i = $len;
      } else {
        for ( ; $b>10; $b--) {
          last if (substr($v,$i+$b,1) =~ /,/) ; # fixme : allowed symbols here
          }
        $b ++; # include last split char
        $w .= substr($v,$i,$b) . "\n$tab";
        $i += $b;
        }
      }
    $v= $w;
    
    $v =~ s/\s*$//;
  }
  return $v;
}

sub togenbank {
  my( $fasta, $cdsaa, $geneGFF, $exonGFF, $cdsGFF)= @_;

  my ($gref, $gsrc, $gtype, $gstart,$gstop,$gscore,$gstrand,$gph,$gattr)= @{ $geneGFF->[0] };

  ## fixme: if possible want expand to EXCLUDE nearby genes, i.e. only intergene regions
  my($offa,$offb)= (@offset>1) ? @offset : (-$GENE_EXPAND, $GENE_EXPAND);
  my($spanstart, $spanstop) = ($gstart + $offa, $gstop + $offb); # add optional expansion
  if($spanstart<1) { $spanstart= 1; }
  my $genedna  = get_dna( $fasta, $gref, $spanstart, $spanstop);
  my $dnalength= length($genedna);
  my $ftvaltab= (" ") x 21;
  my $seqtab= (" ") x 10;
  
  my %gattr= map{ my($k,$v)= split("=",$_,2); $k => $v; } split(";", $gattr);
  
  # assume only one gene/mRNA/CDS entry here ??
  my $mrnaloc= exons2loc($exonGFF,$spanstart-1);
  my $cdsloc = exons2loc($cdsGFF,$spanstart-1);
    $mrnaloc= wrapline($mrnaloc,$ftvaltab);
    $cdsloc= wrapline($cdsloc,$ftvaltab);

  my $expandloc="$gref:$spanstart..$spanstop";
  my $geneloc="$gref:$gstart..$gstop:$gstrand";
  my $geneid=   $gattr{ID};
  
  my $date="01-JAN-2008"; # fixme
  my $genecomm= $gattr{Note} || ""; # /comment="$genecomm"
  my $genename= $gattr{Name} || $gattr{symbol} || $gattr{ID};

  $cdsaa =~ s/\*$//; #?
  my $aasize= length($cdsaa);
  my $definition="DEFINITION  gene-span of $genename at $geneloc";

  $cdsaa =~ s/(.{1,60})/$1\n$ftvaltab/g; $cdsaa =~ s/\s*$//;

  # $genedna =~ s/(.{1,60})/$seqtab$1\n/g; # need spacers, line nums? or not required

  $genedna = lc($genedna); #? is this genbank spec
  $genedna =~ s/(.{1,60})/$1\n/g; # need spacers, line nums? or not required
  my @gdna= split "\n",$genedna;
  my $ib=1;
  foreach my $i (0..$#gdna) {
    my $blen= length($gdna[$i]);
    $gdna[$i] =~ s/(.{1,10})/$1 /g;  $gdna[$i] =~ s/ $//;
    $gdna[$i] = sprintf("%9d ",$ib) . $gdna[$i];
    $ib += $blen;
    }
  $genedna= join("\n",@gdna);
 # augustus etraining whines about wrong length dna/source; not right; needs BASE COUNT line

  return <<"GBEOF";
LOCUS       $geneid        $dnalength bp    DNA             UNA       $date
$definition
COMMENT     Converted from gff+fasta with cdsgff2genbank version n.
FEATURES             Location/Qualifiers                    

     source          1..$dnalength
                     /location="$expandloc"

     mRNA            $mrnaloc
                     /gene="$geneid"
                     /source="$gsrc"
                     /location="$geneloc"
                     /symbol="$genename"

     CDS             $cdsloc
                     /gene="$geneid"
                     /aa_size=$aasize
                     /translation="$cdsaa"

BASE COUNT     0 a   0 c   0 g   0 t
ORIGIN
$genedna
//
GBEOF
}


sub cds2seq {
  my( $fasta, $ref, $cdsA, $seqtype)= @_;
  #old# my($refdna, $cdsA, $seqtype)= @_;
  ## return "" unless(ref $refdna);
  ## assume cdsA are all/only cds exon set for one gene/mrna
  
  my $cstrand= $cdsA->[0]->[6];
  my $isrev= ($cstrand eq '-' || $cstrand < 0);
  
  my @cds= @$cdsA;   # sort by start
  my (@geneGFF, @exonGFF, @cdsGFF);
  if( $seqtype =~ /genbank/) {
    foreach (@cds) {
      if($_->[2] =~ /mRNA|gene/) { push(@geneGFF,$_);  }
      elsif($_->[2] =~ /exon/) { push(@exonGFF,$_);  }
      else {  push(@cdsGFF,$_); } # save, dont sort? dont reverse
      }
    @cds= @cdsGFF;
    }
    
  if ($isrev) { @cds= sort{ $b->[3] <=> $a->[3] } @cds; } # end 1st
  else { @cds= sort{ $a->[3] <=> $b->[3] } @cds; } # start 1st
  my $nt_length= 0;
  my $ispartial= 0;
  my $phase0= 0;
  my $cdsdna= ""; my @exondna=();
  my $id="";
  my ($gref,$gstart,$gstop,$gtype,$gsrc,$gstrand,$gattr)= ("")x10;
  my $header="";
  my $issubset= ($seqtype =~ /start|end/i or $seqtype =~ /intron$/i);
   $issubset=1 if($seqtype =~ /exon/);  # dnaexon ??
 
  foreach my $ix (0 .. $#cds)  {
    my($ref,$src,$type,$start,$stop,$score,$strand,$phase,$attr)= @{$cds[$ix]};
    
    my($rstart,$rstop)= ($start,$stop);
    my($offa,$offb)= @offset;
    
    if($seqtype =~ /intron$/i) { 
      next unless($ix >=0 && $ix < $#cds);
      $src.="-intron"; 
      my($ref1,$src1,$type1,$start1,$stop1)= @{$cds[$ix+1]};
      if($isrev) {
        ($start,$stop)= ($start - 1, $stop1 + 1); # dang is this right?
      } else {
        ($start,$stop)= ($stop + 1, $start1 -1); 
      }
      
    } elsif($issubset && @offset) { 
      if($offa>$offb) { ($offa,$offb)= ($offb,$offa); }
      ($start,$stop)= ($stop,$start)  if($isrev);
      ($offa,$offb) = (-$offa,-$offb) if($isrev);
      if ($seqtype =~ /genestart|intronend/i) { 
        ($start,$stop)= ($start + $offa, $start + $offb); 
        }
      elsif ($seqtype =~ /geneend|intronstart/i) { 
        ($start,$stop)= ($stop + $offa, $stop + $offb); 
        }
      }
    
    #?? fixme for offset err?
    if($start>$stop) { ($start,$stop)= ($stop,$start); }

    if($ix == 0) { 
         if($seqtype =~ /genestart/i) {  $src.="-g5"; }
      elsif($seqtype =~ /geneend/i) {  $src.="-g3"; }
      elsif($seqtype =~ /intronstart/i) { $src.="-in5"; }
      elsif($seqtype =~ /intronend/i) { $src.="-in3"; }
      
      $id=$attr; $id =~ s/^(ID|Parent)=//; $id =~ s/[;,\s].*$//;
      $gattr=$attr; if($gattr=~s/^[^;]+;// && $gattr=~/\w/) {$gattr.=";";} else {$gattr="";}
      ($gref,$gstart,$gstop,$gsrc,$gtype,$gstrand)= ($ref,$start,$stop,$src,$type,$strand);
      $ispartial= 0;
      $phase0= $phase || 0; #? use for atg ?
    }
      
    my($bstart,$blen)= ($start - 1, $stop - $start + 1);

    if($seqtype =~ /genestart/i) { next unless($ix == 0); }
    elsif($seqtype =~ /geneend/i) { next unless($ix == $#cds);  }
    elsif($seqtype =~ /intronstart/i) { next unless($ix >=0 && $ix < $#cds); }
    elsif($seqtype =~ /intronend/i) { next unless($ix >0 && $ix <= $#cds); }
     
    # my $exondna  = substr( $refdna->seq(), $bstart, $blen);
    my $exondna  = get_dna( $fasta, $ref, $start, $stop);
    unless($exondna) {
      warn "# missing $ref:$start-$stop\n" if($verbose);
      return ; # should write header note
      }
      
    if($isrev) {  
      $exondna = reverse $exondna;
      $exondna =~ tr/gatcGATC/ctagCTAG/;
      }

#    if($issubset && $casechange) {
#       # $offa < $offb ... both may be neg or pos.
#       my $ibeg= ($rstart - $bstart+1);
#       my $iin=  ($rstop - $rstart);
#       my $iend= ($stop - $rstop);
#       my $xb= substr($exondna,0,$ibeg);
#       my $xi= substr($exondna,$ibeg,$iin-$ibeg);
#       my $xe= substr($exondna,$iin,$iend-$iin);
#       $xb =~ tr/A-Z/a-z/;
#       $xi =~ tr/a-z/A-Z/;
#       $xe =~ tr/A-Z/a-z/;
#       $exondna = $xb . $xi . $xe;
#      }

    $cdsdna.= $exondna;
    ## ^^ FIXME: don't append like this for $issubset ; need each exon/intron separate fasta entry
    push(@exondna,$exondna);
    
    if( $ix == 0 ) { # only for cds types ...
      ## my $inc5= 0;
      ## 1st exon; find start ATG; ** only need 3 bases at start, not all
      unless($issubset) {
      my $atg= substr($exondna, $phase0, 3);
      $ispartial=1 unless($atg =~ /ATG/i);
      }

    } else {
      $gstart=$start if($start<$gstart);
      $gstop= $stop if($stop>$gstop);
    }
  }
    
  my $minlen= $MIN_CDS_LEN;  
  if($seqtype =~ /aa|prot|genbank/i) { # add trans for $togenbank ?
    my($inc5,$protaa,$aascore);
    if(! $ispartial || ( $trustphase && $phase0 =~ /\d/ )) {
      ($inc5,$protaa,$aascore) = getBestFrame( $cdsdna, $id, ! $ispartial, $phase0); 
      }
    else { 
      ($inc5,$protaa,$aascore) = getBestFrame( $cdsdna, $id, ! $ispartial); 
      }
    $cdsdna= $protaa;
    $minlen= $MIN_AA_LEN;  
    }
  elsif($issubset) { 
    $minlen=1;
    }
    
  my $slen= length($cdsdna);
  return "" if($slen < $minlen);
  
  if($seqtype =~ /genbank/) {
    @cds= sort{ $a->[3] <=> $b->[3] } @cds; # always forward sort for genbank
    @exonGFF= sort{ $a->[3] <=> $b->[3] } @exonGFF; # always forward sort for genbank
    return togenbank($fasta, $cdsdna, \@geneGFF, \@exonGFF, \@cds);
  }
  
  if($issubset) {
    my $i=0; my $fa="";
    foreach my $ex (@exondna) {
      $i++;
      $slen= length($ex);
      #?? want exon/part loc instead of gloc ?
      $header=">$id.$i loc=$gref:$gstart-$gstop:$gstrand;type=$gtype.$gsrc;${gattr}len=$slen";
      # $ex .= "\n"; #$ex =~ s/(.{1,50})/$1\n/g;
      if($prettyfa) { $ex =~ s/(.{1,50})/$1\n/g; } else { $ex .= "\n";  }
      $fa .= $header."\n".$ex;
      }
    return $fa;
  } else {
    $header=">$id loc=$gref:$gstart-$gstop:$gstrand;type=$gtype.$gsrc;${gattr}len=$slen";
    $header .= ";partial_gene=true" if ($ispartial);
    return $header."\n" if($debug);
    if($prettyfa) { $cdsdna =~ s/(.{1,50})/$1\n/g; } else { $cdsdna .= "\n";  }
    return $header."\n".$cdsdna;
  }
}


## translation methods

my @s5CodonTable = ();
BEGIN{
 @s5CodonTable = (
	 [
		 ['K','N','K','N','X',],
		 ['T','T','T','T','T',],
		 ['R','S','R','S','X',],
		 ['I','I','M','I','X',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['Q','H','Q','H','X',],
		 ['P','P','P','P','P',],
		 ['R','R','R','R','R',],
		 ['L','L','L','L','L',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['E','D','E','D','X',],
		 ['A','A','A','A','A',],
		 ['G','G','G','G','G',],
		 ['V','V','V','V','V',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['*','Y','*','Y','X',],
		 ['S','S','S','S','S',],
		 ['*','C','W','C','X',],
		 ['L','F','L','F','X',],
		 ['X','X','X','X','X',],
	],
	 [
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
		 ['X','X','X','X','X',],
	],

);
}

sub ibase {
  my $c= substr($_[0],$_[1],1);
  return 0 if ($c eq 'A');
  return 1 if ($c eq 'C');
  return 2 if ($c eq 'G');
  return 3 if ($c eq 'T');
  return 4;
}  
  
sub translate {
  my($cds, $offset)= @_;
  $cds = uc($cds); ## fix chars ??
  my $aa="";
  my $aa_length = int((length($cds) - $offset) / 3);
	for (my $i = 0; $i < $aa_length; $i++) {
		my $idx = 3 * $i + $offset;
		$aa .= $s5CodonTable[ ibase($cds,$idx)][ ibase($cds,$idx+1) ][ ibase($cds,$idx+2) ];
	}
  return $aa; 
}

sub getBestFrame {
  my($cds, $id, $isfullcds, $usephase)= @_;
  my ($bestscore,$besti,$bestpro)= (-999999,0,"");
  my $ph0= ($usephase =~ /\d/) ? $usephase : 0;
  my $nframe= ($isfullcds ? $ph0+1 : 3);
  for (my $i= $ph0; $i<$nframe; $i++) {
    my $pro= translate( $cds, $i );
    my $score = $pro =~ tr/*/*/; # has_internal_stops($pro);
    #is inner M bad? NO# $score += $pro =~ tr/M/M/;   # has_internal_starts($pro);
    $score *= -3; 
    if (substr($pro,length($pro)-1,1) eq '*') { $score += 4; } # adj internal == end
    if (substr($pro,0,1) eq 'M') { $score += 1; }
    if ($score > $bestscore) { $besti= $i; $bestscore=$score; $bestpro=$pro; }
 warn("# bestFrame[$i,$id]: $score ; $pro \n") if($verbose); # debug  
    last if($score >= 5) # best possible; stop here
    }
  return wantarray ? ($besti,$bestpro,$bestscore) : $besti;
}

__END__

# sub cdsPhase {
#   my($refdna, $cdsA)= @_;
#   ## assume cdsA are all/only cds exon set for one gene/mrna
#   return $cdsA unless(ref $refdna);
#   
#   my $cstrand= $cdsA->[0]->[6];
#   my $isrev= ($cstrand eq '-' || $cstrand < 0);
#   
#   my @cds;   # sort by start
#   if ($isrev) { @cds= sort{ $b->[3] <=> $a->[3] } @$cdsA; } # end 1st
#   else { @cds= sort{ $a->[3] <=> $b->[3] } @$cdsA; } # start 1st
#   my $nt_length= 0;
#   my $ispartial= 0;
#   
#   foreach my $ix (0 .. $#cds)  {
#     my($ref,$src,$type,$start,$stop,$score,$strand,$phase,$attr)= @{$cds[$ix]};
#     my $id=$attr; $id =~ s/^(ID|Parent)=//; $id =~ s/[;,\s].*$//;
#     # do we check exon ordering? use as given in gff?  need 1st .. last, differs for strands
#     
#     if($ix == 0) { 
#       ## 1st exon; find start ATG; ** only need 3 bases at start, not all
#       $ispartial= 0;
#       my($bstart,$blen)= ($start - 1, $stop - $start + 1);
#       #my($bstart,$blen)= ($isrev) ? ($stop-8,8) : ($start-1, 8);
#       my $exondna  = substr( $refdna->seq(), $bstart, $blen);
#       if($isrev) {  
#         $exondna = reverse $exondna;
#         $exondna =~ tr/gatcGATC/ctagCTAG/;
#         }
#       my $inc5= 0;
#       for (; $inc5<=3; $inc5++) {
#         my $atg= substr($exondna, $inc5, 3);
#         last if($atg =~ /atg/i);
#         }
# 
#       ## fixme, if $ispartial probably need check best aa translation frame
#       ## yes; need full cds/all exons and translate() method      
#       if ($inc5 > 2) { 
#         $nt_length = 0;  $ispartial=1; $inc5 = 0; #start not found/incomplete prot ?
#         my $cdsdna= $exondna;
#         foreach my $ex (1 .. $#cds) {
#           my($ref1,$src1,$type1,$start1,$stop1,$score1,$strand1,$phase1,$attr1)= @{$cds[$ex]};
#           my($bstart,$blen)= ($start1 - 1, $stop1 - $start1 + 1);
#           my $exon2dna  = substr( $refdna->seq(), $bstart, $blen);
#           if($strand1 eq '-' || $strand1 < 0) {  
#             $exon2dna = reverse $exon2dna;
#             $exon2dna =~ tr/gatcGATC/ctagCTAG/;
#             }
#           $cdsdna.= $exon2dna;
#           }
#         $inc5 = getBestFrame( $cdsdna, $id);
#         }
#       
#       if ($inc5 == 1) { $nt_length = 2; }
#       elsif ($inc5 == 2) { $nt_length = 1; }
#       else  { $nt_length = 0; }
#     }
#     
#     my($inc5,$inc3,$elength,$frame);
#     $elength = $stop - $start + 1;
# 		$nt_length  += $elength;
# 		$inc3        = $nt_length % 3;
# 		$inc5        = ($elength - $inc3) % 3; # only care about this one
# 		$frame       = ($start + $inc5) % 3;
# 		if ($inc5 == -1) { $inc5 = 2; }
#     
#     my $changed=0;
#     if ($phase eq '.') {  $changed=1; }
#     elsif ($phase ne $inc5 ) { 
#       $changed=2; 
#       warn "# phase change exon[$ix]: $phase => $inc5; $ref:$start-$stop/$strand,$type:$src,$id\n" if $verbose;
#       } 
#     if($changed) { $cds[$ix]->[7]= $inc5; }  
#     if($ispartial && $ix == 0) { $cds[$ix]->[8] .= ";partial_gene=true"; } # 5prime_partial=true; 3prime..
#     }
#     
#   return \@cds;
# }



#item  fasta db

  # dbm based index; not platform independent
  my $dna_db = Bio::DB::Fasta->new('/path/to/fasta/files');

  my $seq     = $db->seq('CHROMOSOME_I',4_000_000 => 4_100_000);
  my $revseq  = $db->seq('CHROMOSOME_I',4_100_000 => 4_000_000);

  my $obj     = $db->get_Seq_by_id('CHROMOSOME_I');
  my $seq     = $obj->seq;
  my $subseq  = $obj->subseq(4_000_000 => 4_100_000);

  #-- or plat-independent (inherits from above) --
  my $lucene = new Bio::DB::GFF::Adaptor::lucene(xxx);
  my $dna_db = new Bio::DB::GFF::Adaptor::LuceneFasta( $fafile, _adaptor => $lucene ); 

#cut



## ATG is universal start codon for euk. nuclear genes

## for start= 0 .. 2, find ATG 
## ? and check aatrans, internal *-stops, aa[0] == 'M', aa[-1]= '*'
## parts from zoeCDS.c of Snap/I.Korf 

	if (exon->strand == '+') {
		c1 = dna->s5[exon->start];
		c2 = dna->s5[exon->start +1];
		c3 = dna->s5[exon->start +2];
		if (c1 == 0 && c2 == 3 && c3 == 2) ef.start = 1; //# ATG
	} else if (exon->strand == '-') {
		c1 = dna->s5[exon->end];
		c2 = dna->s5[exon->end -1];
		c3 = dna->s5[exon->end -2];
		if (c1 == 3 && c2 == 0 && c3 == 1) ef.start = 1; //# TAC (~ATG)
	}
#..........

	best_score = -1000000;
	best_idx   = -1;
	for (i = 0; i <= 2; i++) {
		pro[i] = zoeTranslateDNA(tx->def, tx, i);
		score = - 3 * has_internal_stops(pro[i]);
		if (pro[i]->seq[pro[i]->length -1] == '*') score++;
		if (pro[i]->seq[0] == 'M') score++;
		if (score > best_score) {
			best_score = score;
			best_idx = i;
		}
	}
	
	bt.inc5 = best_idx;
	bt.inc3 = (tx->length - best_idx) % 3;
	bt.aa   = pro[best_idx];
#..........

## find phase == inc5 for all exons:

##	/* label with correct phase and frame */
	if      (cds->inc5 == 1) nt_length = 2;
	else if (cds->inc5 == 2) nt_length = 1;
	else                     nt_length = 0;
		
	elength = 0;
  /*orig* for (i = 0; i < cds->exons->size; i++)  */
  { int i, exstart, exend, exinc;
    if(cds->strand == '-') { exstart= cds->exons->size-1; exend=-1; exinc=-1; }
    else { exstart= 0; exend= cds->exons->size; exinc=1; }
  for (i = exstart; i != exend; i += exinc )
  {
		exon        = cds->exons->elem[i];
		elength     = exon->end - exon->start + 1;
		nt_length  += elength;
		inc3        = nt_length % 3;
		inc5        = (elength - inc3) % 3;
		frame       = (exon->start + inc5) % 3;
		exon->frame = frame;
		exon->inc5  = inc5;
		if (exon->inc5 == -1) exon->inc5 = 2;
		exon->inc3  = inc3;
		if (!zoeVerifyFeature(exon)) {
			zoeWarn("exon does not validate after correcting phase & frame");
			zoeExit("%s", cds->name);
		}
   }
	}
	
