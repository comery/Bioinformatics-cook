use Bio::Tools::GFF;
 
# specify input via -fh or -file
my $gffio = Bio::Tools::GFF->new(-fh => shift, -gff_version => 2);

my $writer = Bio::Tools::GFF->new(-gff_version => 3,
                                  -file        => ">filename.gff3");


my $feature;
# loop over the input stream
while($feature = $gffio->next_feature()) {
    # do something with feature

    #print "$feature\n";
    $writer->write_feature($feature);
}



$gffio->close();
