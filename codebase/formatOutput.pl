#!/usr/bin/perl

#1
$~ = "MYTEXT1";
write;

format MYTEXT1 = 
=================================== 
A *** B *** C *** D *** E *** F ***
=================================== 
.

#2
$~ = "MYTEXT2";
write;
$names = "LILLY";
format MYTEXT2 = 
=================================== 
A *** B *** C *** D *** E *** F ***
$names
=================================== 
.

#3
$~ = "MYTEXT3";
write;

format MYTEXT3 = 
=================================== 
A *** B *** C *** D *** E *** F ***
@<<<<<
$names
=================================== 
.

#4 
open MYFILE, '>', "test.txt" || die $!;  #会在当前目录生成文件test.txt
select(MYFILE);
$~ = "MYTEXT4";
write;
close MYFILE;

format MYTEXT4 = 
=================================== 
A *** B *** C *** D *** E *** F ***
@<<<<<
$names
=================================== 
.


sub write_to_stdout {
  local ($savefile, $saveformat);
  $savefile = select(STDOUT);
  $saveformat = $~;
  $~ = "MYFORMAT";
  write;
  $~ = $saveformat;
  select($savefile);
} 
