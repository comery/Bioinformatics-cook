sub order {
    my $str1 = shift;
	my $lane1 = $1  if ($str1 =~ /(\S+)\/(\d)/);
   #my $order = $lane1 eq $lane2;
 	 return $lane1;
	}
	my $a="FCC4TRMACXX:8:1101:1230:2072#TTAGACAA/1";
	my $b="FCC4TRMACXX:8:1101:1230:2072#TTAGACAA/2";
	my $aa = &order($a);
	my $bb = &order($b);
	print "$aa,$bb";

