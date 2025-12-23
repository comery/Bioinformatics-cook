#!/use/bin/perl

use SVG::TT::Graph::Bar;

my @fields        = qw(Jan Feb Mar);
my @data_sales_02 = qw(12 45 21);

my $graph = SVG::TT::Graph::Bar->new(
  {
        'height' => '500',
		      'width'  => '300',
			        'fields' => \@fields,
					  }
					  );

					  $graph->add_data(
					    {
						      'data'  => \@data_sales_02,
							        'title' => 'Sales 2002',
									  }
									  );

									  open( my $fh, '>', "bar.svg" );
									  select $fh;
									  binmode $fh;
									  print $graph->burn();
									  close($fh);

