#!/usr/bin/perl-w
if (!@ARGV)
   {
    print "Error! No file.\n";
    print "perl circos.pl -fa *.fa -gff *.gff\n -repeat *.RepeatMasker.gff -deep *coverage.txt -outputfile *.png";
    exit;
   }
else
   {
    
      $File_fasta=undef;
      $File_gff=undef;
      $File_deep=undef;
      $File_repeat=undef;
      $outputfile="output.png";
      $value_count=@ARGV;
      for($i=0;$i<$value_count;$i++)
      {  
         if($ARGV[$i]=~/-fa/)
           { 
             $File_fasta=$ARGV[$i+1];
           } 
        
         if($ARGV[$i]=~/-gff/)
           {
             $File_gff=$ARGV[$i+1];
           }
          if($ARGV[$i]=~/-deep/)
           {
             $File_deep=$ARGV[$i+1];
           }
          if($ARGV[$i]=~/-repeat/)
           {
             $File_repeat=$ARGV[$i+1];
           }
           if($ARGV[$i]=~/-outputfile/)
           {
             $outputfile=$ARGV[$i+1];
           }


          $i++;
       }#get input value

##########################creat new gff
       if(!$File_fasta)
          {
              print "Error! No fasta file.\n";
              exit;              
           }
       open NGFF,">>","circos.gff" or die $!;   
       if($File_gff)
       {
          open OGFF,"$File_gff" or die $!;
          while(<OGFF>)
            { 
             print NGFF $_;
             }
        }#exit $File_gff
       close OGFF;
#########################if has fa add CG and gap
        if($File_fasta)
         {
             open FASTA,$File_fasta or die $!;
             <FASTA>;
         if(/>(.*)\s\s\sP/)
           {
         $Seq_id=$1;
         $start=1;
         $end=0;
         $Seq=undef;
         }
          while(<FASTA>)
         {
             chomp;
             if(/>(.*)\s\s\sP/)
               {
                   @array=split(//,$Seq);
                   $length=@array;
                   for($i=0;$i<$length;$i++)
                       {
                               if($array[$i]=~/N/i)
                               {  $start=$i+1;
                                  while($array[$i]=~/N/i)
                                        {   $end=$i+1;
                                            $i++;
                                         }
                                      if($Seq_id)
                                      { print NGFF "$Seq_id\tPATRIC\tgap\t$start\t$end\t.\t+\t0\tCJY_add\n";
                                      }
                                }
                         }
                   $Seq_id=$1;
                   $start=1;
                   $end=0;
                   $Seq=undef;
                   next;
                }
               $Seq.=$_;

          }
         @array=split(//,$Seq);
                   $length=@array;
                   for($i=0;$i<$length;$i++)
                       {  if($array[$i]=~/N/i)
                               {  $start=$i+1;
                                  while($array[$i]=~/N/i)
                                        {   $end=$i+1;
                                            $i++;
                                        }
                                       if($Seq_id)
                                       { print NGFF "$Seq_id\tPATRIC\tgap\t$start\t$end\t.\t+\t0\tCJY_add\n";
                                       }
                                }
                        }
        close FASTA;

####################################add CG to NGFF
     open FASTA,$File_fasta or die $!;
     $line=0;
         while(<FASTA>)
         {    chomp;
             if(/>(.*)\s\s\sP/)
               {
                   if(($line!=0)&&($all!=0))
                   {
                     $GC=$GC_all/$all;
                     print NGFF "$Seq_id\tPATRIC\tGC\t$start\t$end\t.\t+\t$GC\tCJY_add\n";
                    }
                   $Seq_id=$1;
                   $line=0;
                   $GC_all=0;
                   $all=0;
                   $start=1;
                   $end=0;
                   next;

               }
             $line++;
             $countC = $_ =~s/c/c/gi;
             $countG = $_ =~s/g/g/gi;
             $all+=length($_);
             $GC_all += $countC+$countG;
             $end+=length($_);
             if ( $line eq 10)
             {
                  $GC=$GC_all/$all;
                  print NGFF "$Seq_id\tPATRIC\tGC\t$start\t$end\t.\t+\t$GC\tCJY_add\n";
                  $all = 0;
                  $GC_all=0;
                  $line = 0;
                  $start=$end+1;
               }


     }           if($all)
                 { $GC=$GC_all/$all;
                  print NGFF "$Seq_id\tPATRIC\tGC\t$start\t$end\t.\t+\t$GC\tCJY_add\n";
                 }
   close FASTA;


         }#exit $File_fasta
######################### if has deep
         if($File_deep)
         {
          open DEEP,$File_deep or die $!;
          $first_sign=0;
          @Seq=undef;
          while(<DEEP>)
          {
             chomp;
             if(/>(.*)/)
            {
           if($first_sign)
              {
               $size=@Seq;
               $deep=0;
               $count=0;
               for($i=0;$i<$size;$i++)
                   {      $count++;
                       if($count==500)
                        {
                         $avg_deep=$deep/$count;
                         $start=$i-$count+2;
                         $end=$i+1;
                         print NGFF "$seq_id\tPATRIC\tdeep\t$start\t$end\t.\t+\t$avg_deep\tCJYadd\n";
                         $count=0;
                         $deep=0;
                          }
                         $deep+=$Seq[$i];
                    }
                   if($count)
                   {
                     $avg_deep=$deep/$count;
                     $start=$i-$count+1;
                     print NGFF "$seq_id\tPATRIC\tdeep\t$start\t$end\t.\t+\t$avg_deep\tCJYadd\n";
                   }
                  }
               @Seq=undef;
               $first_sign++;
               $seq_id=$1;
           }
         push @Seq,split /\s/;
           }
 ###################print the last one         
              $size=@Seq;
               $deep=0;
               $count=0;
               for($i=0;$i<$size;$i++)
                   {      $count++;
                       if($count==500)
                        {
                         $avg_deep=$deep/$count;
                         $start=$i-$count+2;
                         $end=$i+1;
                         print NGFF "$seq_id\tPATRIC\tdeep\t$start\t$end\t.\t+\t$avg_deep\tCJYadd\n";
                         $count=0;
                         $deep=0;
                          }
                         $deep+=$Seq[$i];
                    }
                   if($count)
                   {
                     $avg_deep=$deep/$count;
                     $start=$i-$count+1;
                     print NGFF "$seq_id\tPATRIC\tdeep\t$start\t$end\t.\t+\t$avg_deep\tCJYadd\n";
                   }
       


         }#exit $File_deep
###############################if has repeat add to NGFF        
         if($File_repeat)
         {
            open REPEAT,$File_repeat or die $!;
         <REPEAT>;
        # <REPEAT>;<REPEAT>;
         while(<REPEAT>)
         {

            my($Seq_id,$start,$end)=(split /\t/)[0,3,4];
            print NGFF "$Seq_id\tPATRIC\trepeat\t$start\t$end\t.\t+\t0\tCJY_add\n";
         #  if(/.*(sid\S*)\s*(\d\S*)\s*(\d\S*)\s*/)
         #   { 
         #     print GFF "$1\tPATRIC\trepeat\t$2\t$3\t.\t+\t0\tCJY_add\n";
         #    }

         }
   close REPEAT;
           
         }#exit $File_repeat
        
 close NGFF;


###################################compelete gff
     open IN,"<$File_fasta" or die " fault:$!\n";
     open OUT,">>","karyotype.txt"or die $!;
    
     $count=0;
     $all_length=0;
     @file_in=undef;
     while(<IN>)
     {
      chomp;
      if(/>(.*)\s\s\sP/)
         {
           if($count)
            {
              push @file_in,"chr - $Seq_id $Seq_id 0 $length blue \n";
             }
           $Seq_id=$1;
           $length=0;
           $count++;
           next;
         }
       else
         {
           $length+=length($_);
           $all_length+=length($_);
         }
      }
      push @file_in,"chr - $Seq_id $Seq_id 0 $length blue \n";
      $left_space=int($all_length/60);
      $right_space=int($all_length/160);
      push @file_in,"chr - left id 0 $left_space white \n";
     # unshift @file_in,"chr - right id 0 $right_space white \n";
      print OUT @file_in;
      close IN;
      close OUT;
####################finish read fa, creat karyotype.txt
###################


###################begin read gff
      open GFF,"circos.gff" or die "fault:$!\n";
      %Filename=undef;
      while(<GFF>)
      {
         chomp;
         if(/^#/)
           {
             next;
            }
         else
            {      my($Seq_id,$type,$start,$end,$content)=(split /\t/)[0,2,3,4,7];
                   if($type=~/.*RNA/)
                     {
                       $type="ncR";
                     }
                   if($type=~/pseudogene/)
                     {
                       next;
                      }
                   if($type=~/repeat/)
                     {
                       $type="RP";
                      }
                    if($type=~/deep/)
                     {
                       $type="DP";
                      }

                   if(!exists $Filename{$type})
                      {  $Filename{$type}=1;
                       }
                   else
                       {
                          $Filename{$type}++;
                       }

                   open OUT,">>","$type.txt"or die $!;

                   if($type=~/GC/)
                     {
                        $tall=$content;
                     }
                   elsif($type=~/DP/)
                       {
                        $tall=$content;

                          }
                   elsif($type=~/CDS/)
                       {
                        $tall=$end-$start+1;
                       }
                     else
                        {
                         $tall=5000;
                        }
                 if($start<=$end)
                  { print OUT "$Seq_id $start $end $tall \n";
                  }
                    close OUT;
                }
           }

      close GFF;
      print "Complete gff convert!\n";


##############################finish gff convert
    open OUT,">>","plot.conf"or die $!;
          print OUT "<plots>\n";
          print OUT "extend_bin=no\n";
          $circle_num=0;
          $File_num=(scalar keys %Filename)-1;
          $space=0.8/$File_num;
          foreach $key (sort keys %Filename )
             {

               if($circle_num==0)
                 {   $circle_num++;
                     next;
                 }

               if($key=~/pseudogene/)
                 {
                     next;
                  }

               open SIGN,">>","$key"."_sign.txt"or die $!;
               print SIGN "left 0 $left_space $key\n";
               close SIGN;

               print OUT "<plot>\n";
               
              
               print OUT "type\t=\thistogram\n";
               print OUT "file\t=\t$key.txt\n";
               $r1=0.97-$space*($circle_num-1);
               $r0=0.97-$space*$circle_num+0.02;
               print OUT "r1\t=\t$r1"."r\n";
               print OUT "r0\t=\t$r0"."r\n";
               if($circle_num<3)
              { print OUT "color\t=\tblack\n";
                print OUT "thickness\t=\t1\n";
               }
                 
               print OUT "fill_under\t=\tyes\n";

              if($key=~/RP/)
               {
                 print OUT "fill_color\t=\tblue\n";
                }
                if($key=~/CDS/)
               {
                 print OUT "fill_color\t=\tred\n";
               }
               if($key=~/GC/)
               {
                  print OUT "fill_color\t=\torange\n";

                }
            #     if($key=~/PG/)
            #   {
            #       print OUT "fill_color\t=\tgreen\n";
            #    }
                   if($key=~/ncR/)
               {
                   print OUT "fill_color\t=\tgreen\n";

                }
                if($key=~/DP/)
               {
                   print OUT "fill_color\t=\tlpurple\n";

                }
               if($key=~/gap/)
               { 
                   print OUT "fill_color\t=\tblack\n";

                }


              # print OUT "background\t=\tyes\n";
              # print OUT "background_color\t=\tvvlgrey\n";
                print OUT "axis\t=\tyes\n";
                print OUT "axis_color\t=\tlgrey\n";
                print OUT "axis_thickness\t=\t2\n";

               if($key=~/GC/)
                  {
                    $axis_spacing=0.2;
                    $min=0;
                    $max=1;
                   }
                elsif($key=~/DP/)
                  {
                    $axis_spacing=40;
                      $min=0;
                    $max=200;
                  }

                 else
                   {
                    $axis_spacing=1000;
                    $min=0;
                    $max=5000;
                    }
               print OUT "min\t=\t$min\n";
               print OUT "max\t=\t$max\n";
               print OUT "axis_spacing\t=\t$axis_spacing\n";
               print OUT "</plot>\n";
               print OUT "\n";
               print OUT "<plot>\n";
               print OUT "type\t=\ttext\n";
               print OUT "file\t=\t$key"."_sign.txt\n";
               print OUT "r1\t=\t$r1"."r\n";
               print OUT "r0\t=\t$r0"."r\n";
               print OUT "color\t=\tblack\n";
               print OUT "label_font\t=\tbold\n";
               print OUT "label_size\t=\t30p\n";
               print OUT "rpadding\t=\t0p\n";
               print OUT "padding\t=\t0p\n";
               print OUT "label_rotate\t=\tno\n";
               print OUT "</plot>\n";
               print OUT "\n";
               $circle_num++;
             }
       print OUT "</plots>\n";
       close OUT;

#########################finish plot.conf

#########################creat circos.conf
          open CONF,">>","circos.conf"or die $!;
          print CONF "karyotype=karyotype.txt\n";
          print CONF "<ideogram>\n<spacing>\ndefault=0.001r\n</spacing>\nradius=0.90r\nthickness=60p\nfill=yes\nstroke_color=white\nstroke_thickness=2p\nshow_label= yes\nlabel_font=default\nlabel_radius=1r+75p\nlabel_size=10p\nlabel_parallel=no\n</ideogram>\n";
          print CONF "show_ticks=yes\nshow_tick_labels=yes\nchromosomes_units = 1000\n<ticks>\nchromosomes=-left;-right\nradius     = 1r\ncolor      = black\nthickness  = 2p\nmultiplier = 1e-3\n<tick>\nspacing    = 1u\nsize       = 3p\ncolor      = lgrey\nshow_label = no\n</tick>\n<tick>\nspacing    = 5u\nsize       = 5p\ncolor      = dgrey\nshow_label = yes\nlabel_size = 8p\nlabel_offset = 0p\nformat       = \%d\n</tick>\n<tick>\nspacing       = 10u\nsize          = 8p\nshow_label    = yes\ncolor         = black\nlabel_size    = 10p\nlabel_offset  = 5p\nformat=\%d\nsuffix        = kb\n</tcik>\n</ticks>\n";
         print CONF "<image>\n<<include //usr/local/Cellar/circos/0.67-7/libexec/etc/image.conf>>\n</image>\n";
         print CONF "<colors>\n<<include //usr/local/Cellar/circos/0.67-7/libexec/etc/colors.conf>>\n</colors>\n";
         print CONF "<fonts>\n<<include //usr/local/Cellar/circos/0.67-7/libexec/etc/fonts.conf>>\n</fonts>\n";
         print CONF "chromosomes_display_default=yes\nchromosomes_scale=/./=1rn\n";
         print CONF "<<include //usr/local/Cellar/circos/0.67-7/libexec/etc/housekeeping.conf>>\n";
         print CONF "<<include plot.conf>>\n";

         close CONF;

system "perl /usr/bin/circos -conf circos.conf -outputfile $outputfile";
#system "rm *.txt";
#system "rm *.conf";
#system "rm *.gff";
   }#correct input value
