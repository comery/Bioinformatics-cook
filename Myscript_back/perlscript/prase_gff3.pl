se strict;
use warnings;
 
use Bio::Tools::GFF;
use Bio::DB::Fasta;
 
my $gff3_file=$ARGV[0];  #gff3 input file
my $genome=$ARGV[1]; #genome sequence directory
my $gffio = Bio::Tools::GFF -> new(-file =>$gff3_file , -gff_version => 3);
my $db    = Bio::DB::Fasta  -> new($genome);
 
my $i=0;
my $first=0;
my $strand=0;
my $seq='';
while(my $feature = $gffio->next_feature()) {
# Sometimes, the gff3 file format is a little different with the standard format,
# So that the keys of the hashes are different.
# Use the following lines of masked code to see the keys of the hashes. Some Change may be needed.
#  $i++;
#  if ($i==2){
#    my @a=keys %{$feature};
#    print "@a\n";
#    my @b=keys %{$feature->{_location}};
#    print "@b\n";
#    print "$feature->{_primary_tag}\n";
#  }
  if($feature->{_primary_tag}=~/mrna/i){
    $first++;
    if($first==1){
      $seq='';
    }else{
      print_sequence($seq,80,$strand);
      $seq='';
    }
    print ">$feature->{_gsf_tag_hash}->{Name}->[0]|$feature->{_gsf_tag_hash}->{ID}->[0]\n";
  }elsif($feature->{_primary_tag}=~/cds/i){
    my $seq_temp=$db->seq($feature->{_gsf_seq_id}, $feature->{_location}->{_start}=>$feature->{_location}->{_end});
    if($feature->{_location}->{_strand}=~/-/){  # to know the strand
      $seq=$seq_temp.$seq;
      $strand=0;
    }else{
      $seq.=$seq_temp;
      $strand=1;
    }
  }
}
print_sequence($seq,80,$strand);
 
$gffio->close();
 
sub print_sequence {
  my($sequence, $length,$strand) = @_;
  if($strand==0){  # if the sequence is on the minus strand, get its reverse-complement counterpart
    $sequence=~tr/ACGTacgt/TGCAtgca/;
    $sequence=reverse $sequence;
  }
  for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ) {
    print substr($sequence, $pos, $length),"\n";
  }
}
