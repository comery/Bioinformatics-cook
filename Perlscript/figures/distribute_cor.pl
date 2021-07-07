#!/usr/bin/perl

=head1 Name

distribute_cor.pl  --  caculate correlation distribution for distribute_svg.pl

=head1 Description

This program is used to prepare data for distribute_svg.pl, to draw correlation
distribution figure for two variables with a given set of number pairs. It can
caculate in two ways: in the same bin size, or in the same data number.  

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2007-2-4

=head1 Usage
 
  % distribute_cor.pl <infile | STDIN>

  --binsize <num>      caculate correlation in specified bin size
  --binnum  <num>      caculate correlation in specified data number
  --minborder <num>    specify the minimum X-axis value cutoff 
  --maxborder <num>    specify the maxiumu X-axis value cutoff

  --header <str>       output a header in format of distribute_svg.pl
  --color  <str>       set a color for the figure
  --mark   <str>       set a mark for the figure

  --verbose   output verbose information to screen  
  --help      output help information to screen  

  .........            all setting items in distribute_svg.pl can be used here

=head1 Exmple

cat chr9.fa.repeat.density | perl ../bin/distribute_cor.pl -binnum 1 -xstart 0 -xend 250000000 -xscale "/1000000" -ystart 0 -yend 1 -yscale "/0.01" -x "chromosome coordinate, Mb" -y "percent of repeat" -note "human chr9" -header point -color blue -pointsize 2  > chr9.fa.repeat.density.lst;

=cut


use strict;
use Getopt::Long;
use Data::Dumper;

my ($BinSize,$BinNum,$MinBorder,$MaxBorder);
my ($Header,$Color,$Mark);
my ($Verbose,$Help);

my ($Width, $Height, $WholeScale);
my ($MarkPos, $MarkScale, $MarkNoBorder, $MarkStyle);
my ($FontSize,$FontFamily);
my ($Note, $X, $Y, $Note2); 
my ($Xstart, $Xstep, $Xend, $XCut, $XScale); 
my ($Ystart, $Ystep, $Yend, $YCut, $YScale); 
my ($XScaleDiv, $YScaleDiv); 
my ($LineWidth); 
my ($PointSize, $Noconnect); 
my ($XUnit, $UnitPer, $MovePer, $OffsetPer);

GetOptions(
	"binsize:f"=>\$BinSize,
	"binnum:n"=>\$BinNum,
	"minborder:f"=>\$MinBorder,
	"maxborder:f"=>\$MaxBorder,

	"header:s"=>\$Header,
	"color:s"=>\$Color,
	"mark:s"=>\$Mark,

	"verbose"=>\$Verbose,
	"help"=>\$Help,
	
	"Width:n"=>\$Width,
	"Height:n"=>\$Height, 
	"WholeScale:f"=>\$WholeScale, 
	
	"MarkPos:s"=>\$MarkPos, 
	"MarkScale:f"=>\$MarkScale, 
	"MarkNoBorder:n"=>\$MarkNoBorder, 
	"MarkStyle:s"=>\$MarkStyle, 
	
	"FontSize:n"=>\$FontSize, 
	"FontFamily:s"=>\$FontFamily, 
	
	"Note:s"=>\$Note,  
	"X:s"=>\$X,  
	"Y:s"=>\$Y,  
	"Note2:s"=>\$Note2, 
	
	"Xstart:f"=>\$Xstart, 
	"Xstep:f"=>\$Xstep, 
	"Xend:f"=>\$Xend, 
	"XCut:n"=>\$XCut,  
	"XScale:s"=>\$XScale, 
	
	"Ystart:f"=>\$Ystart, 
	"Ystep:f"=>\$Ystep, 
	"Yend:f"=>\$Yend, 
	"YCut:n"=>\$YCut,  
	"YScale:s"=>\$YScale,  
	
	"XScaleDiv:n"=>\$XScaleDiv,  
	"YScaleDiv:n"=>\$YScaleDiv,  
	

	"LineWidth:n"=>\$LineWidth, 

	"PointSize:n"=>\$PointSize,  
	"Noconnect:n"=>\$Noconnect,  

	"XUnit:f"=>\$XUnit, 
	"UnitPer:f"=>\$UnitPer,    
	"MovePer:f"=>\$MovePer,  
	"OffsetPer:f"=>\$OffsetPer

);


die `pod2text $0` if ($Help);

my @data;
my $total;
my %Xpre;
my %X;
my $output;

while (<>) {
	push @data,[$1,$2] if(/([-\d\.eE]+)\s+([-\d\.eE]+)/);
}
@data = sort {$a->[0]<=>$b->[0]} @data;
$MinBorder = $data[0][0] unless(defined $MinBorder);
$MaxBorder = $data[-1][0] unless(defined $MaxBorder);

print STDERR "read data done\n" if($Verbose);

if (defined $BinSize) {
	
	$BinSize = ($MaxBorder - $MinBorder) / 50 unless($BinSize);

	##skip numbers lower than $MinBorder
	my $data_pos = 0;
	foreach my $p (@data) {
		if ($p->[0] < $MinBorder) {
			$data_pos++;
		}else{
			last;
		}
	}

	my ($bin_start,$bin_end,$bin_mid);
	for ($bin_start=$MinBorder; $bin_start<$MaxBorder; $bin_start+=$BinSize) {
		$bin_end = $bin_start + $BinSize;
		$bin_mid = $bin_start + $BinSize/2;
		while ($data_pos<scalar(@data)) {
			last if($data[$data_pos][0] >= $bin_end);
			push @{$Xpre{$bin_mid}},$data[$data_pos];	
			$total++;
			$data_pos++;
		}
	}

	## include numbers equal $MaxBorder
	while ($data_pos<scalar(@data)) {
		last if($data[$data_pos][0] > $MaxBorder);
		push @{$Xpre{$bin_mid}},$data[$data_pos];
		$total++;
		$data_pos++;
	}

	foreach my $bin_mid (sort {$a<=>$b} keys %Xpre) {
		my $bin_p = $Xpre{$bin_mid};
		my ($unit_x,$unit_y,$unit_num);
		foreach my $p (@$bin_p) {
			$unit_x += $p->[0];
			$unit_y  += $p->[1];
			$unit_num++;
		}
		$X{$bin_mid} = $unit_y / $unit_num;
	}

}


if (defined $BinNum) {
	
	$BinNum = 1 unless($BinNum);
	
	my ($unit_x,$unit_y,$unit_num);
	foreach my $p (@data) {
		next if($p->[0] < $MinBorder);
		last if($p->[0] > $MaxBorder);
		$unit_x += $p->[0];
		$unit_y += $p->[1];
		$unit_num ++ ;
		if ($unit_num == $BinNum) {
			$X{$unit_x/$unit_num} = $unit_y/$unit_num;
			$unit_x = 0;
			$unit_y = 0;
			$unit_num = 0;
		}
	}
}


print STDERR "bin caculate done\n" if($Verbose);


if (defined $Header) {
	&header_default();
	$output .= &header_output();
}

$output .= "\nColor: $Color\n" if(defined $Color);
$output .= "Mark: $Mark\n" if(defined $Mark);
foreach my $bin_mid (sort {$a<=>$b} keys %X) {
	$output .= "$bin_mid: $X{$bin_mid}\n";
}

print $output;



####################################################
################### Sub Routines ###################
####################################################



sub header_default{
	
	$Width = 640 unless(defined $Width);
	$Height = 480 unless(defined $Height);
	$WholeScale = 0.8 unless(defined $WholeScale);
	
	$MarkPos = 'tr' unless(defined $MarkPos);
	$MarkScale = 0.8 unless(defined $MarkScale);
	$MarkNoBorder = 0  unless(defined $MarkNoBorder);
	$MarkStyle = 'v' unless(defined $MarkStyle);
	
	$FontSize = 46 unless(defined $FontSize);
	$FontFamily = 'ArialNarrow-Bold'  unless(defined $FontFamily);
	
	unless(defined $Note){
		$Note = 'Correlation';
	}
	$X = 'Variable 1' unless(defined $X);
	$Y = 'Variable 2' unless(defined $Y);
	$Note2 = '' unless(defined $Note2);
	
	$Xstart = $MinBorder unless(defined $Xstart);
	$Xend = $MaxBorder unless(defined $Xend);
	$Xstep = ($Xend - $Xstart) / 5 unless(defined $Xstep);
	$XCut = 0 unless(defined $XCut);

	my @temp;
	if(! defined $XScale){	
		$XScale = '';
		for (my $i=$Xstart; $i<=$Xend; $i+=$Xstep) {
			push @temp,$i;
		}
		
	}elsif($XScale =~ /\/([-\d\.Ee]+)/){ ## 除以某个尺度
		$XScale = "";
		for (my $i=$Xstart; $i<=$Xend; $i+=$Xstep) {
			push @temp, $i / $1;
		}

	}else{
		$XScale =~ s/^\s+//;
		$XScale =~ s/\s+$//;
		@temp = split /\s+/,$XScale;
		$XScale = "";
		foreach  (@temp) {
			s/^_/-/;
		}
	}
	
	foreach  (@temp) {
		$XScale .= $_."\n";
	}
	
	unless(defined $Ystart){
		foreach my $bin_mid (sort {$a<=>$b} keys %X) {
			$Ystart = $X{$bin_mid} if($X{$bin_mid} < $Ystart || !$Ystart);
		}
	}

	unless(defined $Yend){
		foreach my $bin_mid (sort {$a<=>$b} keys %X) {
			$Yend = $X{$bin_mid} if($X{$bin_mid} > $Yend || !$Yend);
		}
	}
	$Ystep = ($Yend - $Ystart) / 5 unless(defined $Ystep);
	$YCut = 0 unless(defined $YCut);
	
	my @temp;
	if(! defined $YScale){	
		
		$YScale = '';
		for (my $i=$Ystart; $i<=$Yend; $i+=$Ystep) {
			push @temp,$i;
		}
		
	}elsif($YScale =~ /\/([-\d\.Ee]+)/){ ## 除以某个尺度
		$YScale = "";
		for (my $i=$Ystart; $i<=$Yend; $i+=$Ystep) {
			push @temp, $i / $1;
		}

	}else{
		$YScale =~ s/^\s+//;
		$YScale =~ s/\s+$//;
		@temp = split /\s+/,$YScale;
		$YScale = "";
		foreach  (@temp) {
			s/^_/-/;
		}
	}
	foreach  (@temp) {
		$YScale .= $_."\n";
	}

	$XScaleDiv = 1 unless(defined $XScaleDiv);
	$YScaleDiv = 1 unless(defined $YScaleDiv);
	

	$LineWidth = 3 unless(defined $LineWidth);

	$PointSize = 3 unless(defined $PointSize);
	$Noconnect = 1 unless(defined $Noconnect);

	$XUnit = 1  unless(defined $XUnit);
	$UnitPer = $BinSize  unless(defined $UnitPer);
	$MovePer = -$BinSize/2 unless(defined $MovePer);
	$OffsetPer = 0 unless(defined $OffsetPer);

}


sub header_output{
	my ($common,$line,$point,$rect);

	$common = <<HEADER;
Width:$Width
Height:$Height
WholeScale:$WholeScale
MarkPos:$MarkPos
MarkScale:$MarkScale
MarkNoBorder:$MarkNoBorder
MarkStyle:$MarkStyle
FontSize:$FontSize
FontFamily:$FontFamily
Note:$Note 
X:$X
Y:$Y
Xstart:$Xstart
Xstep:$Xstep 
Xend:$Xend
XCut:$XCut
XScale:
$XScale:End
Ystart:$Ystart
Ystep:$Ystep 
Yend:$Yend
YCut:$YCut
YScale:
$YScale:End
XScaleDiv:$XScaleDiv
YScaleDiv:$YScaleDiv
Note2:
$Note2:End
HEADER
	
	
	$line = <<HEADER;
Type:Line
LineWidth:$LineWidth
$common:End

HEADER

	$point =  <<HEADER;
Type:Point
PointSize:$PointSize
Noconnect:$Noconnect
$common:End

HEADER

	$rect =  <<HEADER;
Type:Rect
XUnit:$XUnit
UnitPer:$UnitPer     
MovePer:$MovePer  
OffsetPer:$OffsetPer
$common:End

HEADER

	return $line if( $Header =~ /line/i);
	return $point if( $Header =~ /point/i);
	return $rect if( $Header =~ /rect/i);
}



