#!/usr/bin/perl

open (STDOUT, ">file1");     #写
open (STDERR, ">&STDOUT");   #流重写向
open (READWRITE, "+>file1"); #既或读又或写

$| = 1;   #文件缓冲,如果$|为非零值则不使用缓冲
select (STDERR); #当未调用select函数时，$|影响当前缺省文件
if (eof)    #文件结束，通常可在while中使用

while ($line = <>)  #按行读取, <>中可以写句柄，默认STDIN




