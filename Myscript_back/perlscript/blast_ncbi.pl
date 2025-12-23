use strict;
use warnings;
die "Usage:\n\tperl $0 query.fa outdir" unless (@ARGV ==2);
use Bio::Tools::Run::RemoteBlast;
my $prog='blastn';
my $db='nt';
my $e_value='1e-10';
my $factory = Bio::Tools::Run::RemoteBlast->new('-prog' => $prog,'-data' => $db,'-expect' => $e_value,'-readmethod' => 'SearchIO' );
my $v = 1;
my $str = Bio::SeqIO->new(-file=>shift , -format => 'fasta' );
my $outdir = shift;
mkdir $outdir unless (-d $outdir);
open LOG, ">$outdir/log";
while (my $input = $str->next_seq())
{
my $r = $factory->submit_blast($input);
print LOG "waiting..." if( $v > 0 );
    while ( my @rids = $factory->each_rid )
{
    foreach my $rid ( @rids ) 
  {
    my $rc = $factory->retrieve_blast($rid);
     if( !ref($rc) ) 
   {
              if( $rc < 0 ) 
      {
                $factory->remove_rid($rid);
       }
              print LOG "." if ( $v > 0 );
              sleep 1;
            } else 
   {
              my $result = $rc->next_result();
              #save the output
              my $filename = $result->query_name().".out";
			  $filename = $outdir."/"."$filename";
              $factory->save_output($filename);
              $factory->remove_rid($rid);
              print "nQuery Name: ", $result->query_name(), "n";
              while ( my $hit = $result->next_hit )
   {
                next unless ( $v > 0);
                print "thit name is ", $hit->name, "n";
                while( my $hsp = $hit->next_hsp ) 
    {
                  print "ttscore is ", $hsp->score, "n";
                }
              }
            }
          }
        }
}


#-----------------------------------------------------------------#
#this perl script is for throw your query to ncbi database 
#for serching #blast result !
#-----------------------------------------------------------------#
