sub cut2fmt {
	my(@positions) = @_;
	my $template = '';
  	my $lastpos = 1;
  	foreach $place (@positions) {
	$template .= "A" . ($place - $lastpos) . " ";
	$lastpos = $place;
   	}
	$template .= "A*";
return $template;
}
$fmt = cut2fmt(8, 14, 20, 26, 30);
print "$fmtn";
# A7 A6 A6 A6 A4 A*
#-----------------------------
