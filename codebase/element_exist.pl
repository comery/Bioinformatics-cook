#! /usr/bin/perl
use Digest::MD5 qw(md5 md5_hex md5_base64);
my $size = 99999; 
my @array = map MD5->hash($_), 1 .. $size; 
my %hash; @hash{@array} = (1) x @array; 
foreach my $seed (0, int($size/2), 1) 
{ 
     print '-' x 72, "\n\$key = MD5('$seed')\n"; 
     my $key = MD5->hash($seed); 
     cmpthese(-1, { 
            'hash'  => sub { exists $hash{$key} }, 
            'hash+' => sub { @hash{@array} = (1) x @array; exists $hash{$key} }, 
            #'~~'    => sub { $key ~~ @array }, 
            'grep'  => sub { grep $_ eq $key, @array}, 
            'for'   => sub { foreach (@array) { last if $_ eq $key} }, 
     }); 
 } 
 $key = MD5('49999')
