#!usr/bin/perl -w
use strict;
foreach my $file (glob "^..*") {
	my $filename =$file;
	my @neme=split (/\./,$filename);
	my $namestr=$name[1];
	my $newfile=~ s/$filename/RFUNaqnTAARAAPEI-41.$namestr/;
	if (-e $newfile) {
		warn "can't rename $file to $newfile:$newfile exists\n";
		}elsif (rename $file,$newfile) {
		}else {
		warn "rename $file to $newfile failed:$!\n";
		}
	}
