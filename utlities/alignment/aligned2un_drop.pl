open IN,shift;
while(my $id = <IN>) {
	chomp $id;
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/-//g;
	if (length($seq) > 10){
		print "$id\n$seq\n";
	}
}
