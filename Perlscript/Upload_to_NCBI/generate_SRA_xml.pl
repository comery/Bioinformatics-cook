#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use POSIX qw(strftime);

my $version = "V1.5.2 / 140103";
my $XSD_DIR = "$Bin/xsd/SRA_1-5";

my ($out_dir, $validate_xml, $help);
GetOptions( 
	"o=s" => \$out_dir,
	"h"   => \$help);

my $usage=<<USAGE;
Program: generate_SRA_xml.pl (A program generate xml for submit sequencing data to NCBI SRA database)
Version: $version
Contact: Chen Chao <chenchao\@genomics.cn>\n
Usage:   perl generate_SRA_xml.pl [options] <confie_file>

Options: -o  STR   output file directory your want [./]
         -h        print this help information\n
Example: perl generate_SRA_xml.pl -o /share/backup/chenchao/selfbin/SRA_xml configure.txt\n
USAGE
die $usage if (@ARGV < 1 or $help);

$out_dir ||= "./";
my ($config_file) = shift;
open CONFIG, $config_file or die $!;

my ($sra_id, $sra_title, $study_type, $study_abstract, $center_name ); ##study_description
my ($species, $taxon_id, $instrument_model, $base_caller, $sequence_space, $library_construction_protocol, $qtype, $multiplier, $number_of_levels, $quality_scorer); ##expriment_description
my ($fastq_list, $data_mark, $bioproject_id, $library_strategy, $library_source, $library_selection, $design_description); ##difference_option
my ($submitter_name, $submitter_email); ##submission_description
my (%difference_option, %fq_count, %sample_count);
my ($string, $fq_list);

########## deal with the config file

while (<CONFIG>) {
	next if (/^#/ or /^\s/);
	chomp;
	my @F = split /=/;
	die "ERROR: $F[0] is empty maybe have blank space between element and value!!\t" unless $F[1]=~/\w/;
	my $attribute = (split /#/,$F[1])[0];
	$attribute =~ s/=|^\s+|\s+$//;
	

	##study_description
	if ( /^SRA_ID/ ) { $sra_id = $attribute; }
	elsif ( /^SRA_TITLE/ ) { $sra_title = $attribute; }
	elsif ( /^STUDY_TYPE/ ) { $study_type = $attribute; }
	elsif ( /^STUDY_ABSTRACT/ ) { $study_abstract = $attribute; }
	elsif ( /^CENTER_NAME/ ) { $center_name = $attribute; }
	#elsif ( /^CENTER_PROJECT_NAME/ ) { $center_project_name = $F[1]; }

	##expriment_description
	elsif ( /^Species/ ) { $species = $attribute; }
	elsif ( /^TAXON_ID/ ) { $taxon_id = $attribute; }
	elsif ( /^instrument_model/ ) { $instrument_model = $attribute; }
	elsif ( /^BASE_CALLER/ ) { $base_caller = $attribute; }
	elsif ( /^SEQUENCE_SPACE/ ) { $sequence_space = $attribute; }
	elsif ( /^LIBRARY_CONSTRUCTION_PROTOCOL/ ) { $library_construction_protocol = $attribute; }
	elsif ( /^QUALITY_TYPE/ ) { $qtype = $attribute; }
	elsif ( /^MULTIPLIER/ ) { $multiplier = $attribute; }
	elsif ( /^NUMBER_OF_LEVELS/ ) { $number_of_levels = $attribute; }
	elsif ( /^QUALITY_SCORER/ )	{ $quality_scorer = $attribute; }
	##submission_description
	elsif ( /^Submitter_Name/ ) { $submitter_name = $attribute; }
	elsif ( /^Submitter_Email/ ) { $submitter_email = $attribute; }

	##difference_option
	elsif ( /^Fastq_list/ ) {
		$fq_count{$F[0]}++;
		if ( $fq_count{$F[0]} >1 ) {
			$string = "$data_mark\t$bioproject_id\t$library_strategy\t$library_source\t$library_selection\t$design_description";
			$difference_option{$fastq_list} = $string;
			}
		$fastq_list = $attribute;
		}
	elsif ( /^Data_Mark/ ) { $data_mark = $attribute; }
	elsif ( /^BioProject_ID/ ) { $bioproject_id = $attribute; }
	elsif ( /^LIBRARY_STRATEGY/ ) { $library_strategy = $attribute; }
	elsif ( /^LIBRARY_SOURCE/ ) { $library_source = $attribute; }
	elsif ( /^LIBRARY_SELECTION/ ) { $library_selection =  $attribute; }
	elsif ( /^DESIGN_DESCRIPTION/ ) { $design_description = $attribute; }
	}

$string = "$data_mark\t$bioproject_id\t$library_strategy\t$library_source\t$library_selection\t$design_description";
$difference_option{$fastq_list} = $string;

print STDERR "Creating xml code.....\n";
########### study xml
my $studyxml =<<STUDY;
<?xml version="1.0" encoding="UTF-8"?>
<STUDY_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <STUDY alias="$sra_title">
    <DESCRIPTOR>
      <STUDY_TITLE>$sra_title</STUDY_TITLE>
      <STUDY_TYPE existing_study_type="$study_type"></STUDY_TYPE>
      <STUDY_ABSTRACT>$study_abstract</STUDY_ABSTRACT>
      <CENTER_NAME>$center_name</CENTER_NAME>
      <CENTER_PROJECT_NAME>$sra_title</CENTER_PROJECT_NAME>
    </DESCRIPTOR>
  </STUDY>
</STUDY_SET>
STUDY

########### submission_xml
my $current_date = strftime("%Y-%m-%d", localtime(time)); ##get current date 2013-01-18
my $submission_date = "$current_date"."T08:00:00";
my $submission_comment = "submission prepared by $submitter_email";
my $checksum_method = "MD5";

my $submissionxml =<<SUBMISSION;
<?xml version="1.0" encoding="UTF-8"?>
<SUBMISSION_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <SUBMISSION accession="$sra_id" alias="$sra_title" center_name="$center_name" lab_name="$center_name" submission_comment="$submission_comment" submission_date="$submission_date">
    <CONTACTS>
      <CONTACT name="$submitter_name" inform_on_error="$submitter_email" inform_on_status="$submitter_email"/>
    </CONTACTS>
    <ACTIONS>
      <ACTION>
        <ADD schema="experiment" source="$sra_id.experiment.xml"/>
      </ACTION>
      <ACTION>
        <ADD schema="run" source="$sra_id.run.xml"/>
      </ACTION>
    </ACTIONS>
  </SUBMISSION>
</SUBMISSION_SET>
SUBMISSION
##accession --> SRA061649; center_name -- > BGI; handle --> .; lab_name --> BGI; ......
##v1.5 -> delete handle=\".\"; delete submission_id=\"$study_id\";
##$submission .= "\t</FILES>\n\n"; #v1.5 delete
##$submission .= "\t<FILES>\n"; #v1.5 delete
##v1.5 -> delete notes\=\"sample descriptor\"; delete notes\=\"experiment descriptor\"; delete notes=\"runs wrapper\";

########### deal with the Fastq list
my ($subsample,$subexperiment,$subrun);
foreach $fastq_list ( keys %difference_option ) {
	my %count;
	($data_mark,$bioproject_id,$library_strategy,$library_source,$library_selection,$design_description) = split /\t/,$difference_option{$fastq_list};

	open FASTQLIST,$fastq_list or die $!;
	while (<FASTQLIST>) {
		chomp;
		next if (/^#sampleID/ or /^sampleID/ or /^#/);
		my @F = split /\t/;
		for ( my $i=0; $i<@F; $i++ ) { $F[$i]=~s/\s//g } #clean blank space

		my @fastq;
		@fastq = <$F[10]/*fq.gz> or @fastq = <$F[10]*fq>;
		next unless @fastq;
		my ($read1_len,$read2_len) = (0,0);
		if ( $F[6] =~ /\// ) { #pair-end sequence
			($read1_len,$read2_len) = split /\//,$F[6];
			$fq_list .= "$F[12]\t$F[0]\t$F[7]\t$fastq[0]\n$F[13]\t$F[0]\t$F[7]\t$fastq[1]\n";
			}
		else { #single end sequence
			$read1_len = $F[6];
			$fq_list .= "$F[12]\t$F[0]\t$F[7]\t$fastq[0]\n";
			}

		############ sample xml
		my $sample = "$F[0]";
		$sample_count{$sample}++; #sample count
		if ( $sample_count{$sample} == 1 ) {
			$subsample .= &sample_define($sample,$species,$taxon_id) ."\n";
			}

		############ split lane	
		my ($run_date,$instrument_name,$sector) = (split /_/,$F[7])[0,2,3];

		############ experiment xml
		my $library_name = $F[3];
		my $library_mark = "$run_date-$library_name";
		$count{$library_mark}++; #lib count
		my $experiment = "$sample-$data_mark-$library_mark";

		if ( $count{$library_mark} == 1 ) {
			my $nominal_length = $F[4]; #InsertSize

			my $nominal_sdev ;
			$F[5] =~ /-(\d+)\/\+(\d+)/;
			if ( $1 && $2 ) { $nominal_sdev = ($1+$2)/2 } #SDEV
			else { $nominal_sdev = 0 }

			my $orientation = "5'-3'-3'-5'";
			
			my $number_of_reads_per_spot;
			if ( $F[6] =~ /\// ) { $number_of_reads_per_spot= 2 }
			else { $number_of_reads_per_spot= 1 }
			
			my $read_index = "0";
			my $read_class = "Application Read";
			my $read_type = "Forward";
			my $base_coord = "1";
			my $cycle_count = $read1_len + $read2_len;

			$subexperiment .= "  <EXPERIMENT alias=\"$experiment\">\n"; ##alias --> experiment_id
			$subexperiment .= "    <TITLE>$sra_title</TITLE>\n"; ##title --> Multi_sequencing of 99 bladder cancer
			$subexperiment .= "    <STUDY_REF accession=\"$bioproject_id\"/>\n"; ##accession & refname --> SRA061649 & bladder cancer exomes sequencing
			$subexperiment .= "    <DESIGN>\n";
			$subexperiment .= "      <DESIGN_DESCRIPTION>$design_description</DESIGN_DESCRIPTION>\n"; ##DESIGN_DESCRIPTION --> We sequenced the exomes of 99 TCC patients
			$subexperiment .= "      <SAMPLE_DESCRIPTOR accession=\"$sample\"/>\n"; ##refname2 --> $F[7]
			$subexperiment .= "      <LIBRARY_DESCRIPTOR>\n";
			$subexperiment .= "        <LIBRARY_NAME>$library_name</LIBRARY_NAME>\n"; ##LIBRARY_NAME --> $F[5]
			$subexperiment .= "        <LIBRARY_STRATEGY>$library_strategy</LIBRARY_STRATEGY>\n"; ##LIBRARY_STRATEGY --> WGS
			$subexperiment .= "        <LIBRARY_SOURCE>$library_source</LIBRARY_SOURCE>\n"; ##LIBRARY_SOURCE --> GENOMIC
			$subexperiment .= "        <LIBRARY_SELECTION>$library_selection</LIBRARY_SELECTION>\n"; ##LIBRARY_SELECTION --> other
			$subexperiment .= "        <LIBRARY_LAYOUT>\n";
			
			if ( $F[6] =~ /\// ) { #pair-end sequence
				$subexperiment .= "          <PAIRED NOMINAL_LENGTH=\"$nominal_length\" NOMINAL_SDEV=\"$nominal_sdev\"/>\n";
				} ##NOMINAL_LENGTH & NOMINAL_SDEV &ORIENTATION --> 200 & & 5'-3'-3'-5'
				## v1.5 -> delete ORIENTATION=\"$orientation\"
			else { $subexperiment .= "          <SINGLE/>\n"; }
			
			$subexperiment .= "        </LIBRARY_LAYOUT>\n";
			$subexperiment .= "        <LIBRARY_CONSTRUCTION_PROTOCOL>$library_construction_protocol</LIBRARY_CONSTRUCTION_PROTOCOL>\n"; ##LIBRARY_CONSTRUCTION_PROTOCOL --> Standard Solexa protocol
			$subexperiment .= "      </LIBRARY_DESCRIPTOR>\n";
			$subexperiment .= "      <SPOT_DESCRIPTOR>\n";
			$subexperiment .= "        <SPOT_DECODE_SPEC>\n";
			##$subexperiment .= "\t        <NUMBER_OF_READS_PER_SPOT>$number_of_reads_per_spot<\/NUMBER_OF_READS_PER_SPOT>\n"; ##NUMBER_OF_READS_PER_SPOT PE测序需要填两个
			##v1.5 delete
			$subexperiment .= "          <SPOT_LENGTH>$cycle_count</SPOT_LENGTH>\n";
			$subexperiment .= "          <READ_SPEC>\n";
			$subexperiment .= "            <READ_INDEX>$read_index</READ_INDEX>\n"; ##READ_INDEX --> 0
			$subexperiment .= "            <READ_CLASS>$read_class</READ_CLASS>\n"; ##READ_CLASS
			$subexperiment .= "            <READ_TYPE>$read_type</READ_TYPE>\n"; ##READ_TYPE
			$subexperiment .= "            <BASE_COORD>$base_coord</BASE_COORD>\n"; ##BASE_COORD
			##v1.5 delete
			$subexperiment .= "          <\/READ_SPEC>\n";

			if ( $F[6] =~ /\// ) { #pair-end sequence
				$read_index = "1";
				$read_type = "Reverse";
				$base_coord = $read1_len+1;

				$subexperiment .= "          <READ_SPEC>\n";
				$subexperiment .= "            <READ_INDEX>$read_index</READ_INDEX>\n"; ##READ_INDEX --> 1
				$subexperiment .= "            <READ_CLASS>$read_class</READ_CLASS>\n"; ##READ_CLASS
				$subexperiment .= "            <READ_TYPE>$read_type</READ_TYPE>\n"; ##READ_TYPE
				$subexperiment .= "            <BASE_COORD>$base_coord</BASE_COORD>\n"; ##BASE_COORD
				$subexperiment .= "          </READ_SPEC>\n";
				}

			$subexperiment .= "        </SPOT_DECODE_SPEC>\n";
			$subexperiment .= "      </SPOT_DESCRIPTOR>\n";
			$subexperiment .= "    </DESIGN>\n";
			$subexperiment .= "    <PLATFORM>\n";
			$subexperiment .= "      <ILLUMINA>\n";
			$subexperiment .= "        <INSTRUMENT_MODEL>$instrument_model</INSTRUMENT_MODEL>\n"; ##INSTRUMENT_MODEL --> Illumina Genome Analyzer II
			##$subexperiment .= "\t      <CYCLE_SEQUENCE></CYCLE_SEQUENCE>\n";
			##$subexperiment .= "\t      <CYCLE_COUNT>$cycle_count</CYCLE_COUNT>\n"; ##CYCLE_COUNT --> $F[10]*2
			$subexperiment .= "      </ILLUMINA>\n";
			$subexperiment .= "    </PLATFORM>\n";
			$subexperiment .= "    <PROCESSING>\n";
			##$subexperiment .= "\t    <BASE_CALLS>\n";
			##$subexperiment .= "\t      <BASE_CALLER>$base_caller</BASE_CALLER>\n"; ##BASE_CALLER --> Solexa primary analysis
			##$subexperiment .= "\t      <SEQUENCE_SPACE>$sequence_space</SEQUENCE_SPACE>\n"; ##SEQUENCE_SPACE -->Base Space
			##$subexperiment .= "\t    <\/BASE_CALLS>\n";
			##$subexperiment .= "\t    <QUALITY_SCORES qtype=\"$qtype\">\n"; ##qtype  --> other
			###$subexperiment .= "\t      <MULTIPLIER>$multiplier</MULTIPLIER>\n"; ##MULTIPLIER --> 1
			##$subexperiment .= "\t      <NUMBER_OF_LEVELS>$number_of_levels</NUMBER_OF_LEVELS>\n"; ##NUMBER_OF_LEVELS --> 80
			##$subexperiment .= "\t      <QUALITY_SCORER>$quality_scorer</QUALITY_SCORER>\n"; ##QUALITY_SCORER  --> Solexa primary analysis
			##$subexperiment .= "\t    </QUALITY_SCORES>\n";
			$subexperiment .= "    </PROCESSING>\n";
			$subexperiment .= "    <EXPERIMENT_ATTRIBUTES>\n";
			$subexperiment .= "      <EXPERIMENT_ATTRIBUTE>\n";
			$subexperiment .= "        <TAG>center_name</TAG>\n"; ##TAG --> center_name
			$subexperiment .= "        <VALUE>$center_name</VALUE>\n"; ##VALUE --> BGI
			$subexperiment .= "      </EXPERIMENT_ATTRIBUTE>\n";
			$subexperiment .= "    </EXPERIMENT_ATTRIBUTES>\n";
			$subexperiment .= "  </EXPERIMENT>\n\n";
			}

		########### submission_xml
		my ($checksum1, $fastq1, $checksum2, $fastq2);
		if ( $F[6] =~ /\// ) { #pair-end sequence
			$checksum1 = $F[12];
			$fastq1 = basename $fastq[0];
			#$submission .= "\t  <FILE checksum=\"$checksum1\" checksum_method=\"$checksum_method\" filename=\"$fastq1\"/>\n"; #v.15 delete

			$checksum2 = $F[13];
			$fastq2 = basename $fastq[1];
			#$submission .= "\t  <FILE checksum=\"$checksum2\" checksum_method=\"$checksum_method\" filename=\"$fastq2\"/>\n"; #v.15 delete
			}
		else{
			$checksum1 = $F[12];
			$fastq1 = basename $fastq[0];
			#$submission .= "\t  <FILE checksum=\"$checksum1\" checksum_method=\"$checksum_method\" filename=\"$fastq1\"/>\n"; #v.15 delete
			}

		########### run_xml
		my $run = "$sample-$data_mark-$F[7]";
		$run_date =~ s/(\d{2})/$1-/g;
		$run_date =~ s/-$//;
		$run_date = "20$run_date"."T08:00:00";
		my $total_data_blocks = "1";
		my $format_code = "1";
		my $number_channels = "8";
		my $region = "0";

		$sector =~ s/L//;
		my $total_reads = $F[11];
		my $total_spots =$F[11];
		my $filetype = "fastq";
		my $total_bases;
		if ( $F[6] =~ /\// ) { $total_bases = ($read1_len + $read2_len) * $F[11]; }
		else { $total_bases = $read1_len * $F[11]; }

		$subrun .= "  <RUN alias=\"$run\" run_center=\"$center_name\" run_date=\"$run_date\" >\n";
		##v1.5 -> delete instrument_model\=\"$instrument_model\" instrument_name\=\"$instrument_name\" total_data_blocks\=\"$total_data_blocks\"
		##alias --> SRA061649; instrument_model --> Illumina Genome Analyzer II; instrument_name --> $F[3]/FC80VH3ABXX; run_center --> BGI
		$subrun .= "    <EXPERIMENT_REF  refname=\"$experiment\"/>\n"; #refname --> BGI-$F[7]
		$subrun .= "    <DATA_BLOCK>\n";
		##v1.5 -> delete format_code\=\"$format_code\" name\=\"$run\" number_channels\=\"$number_channels\" region\=\"$region\" sector=\"$sector\" total_reads\=\"$total_reads\" total_spots=\"$total_spots\"
		##format_code --> 1; name --> BGI-$F[3]-$F[4]; ......; $total_spots --> $F[11];
		$subrun .= "      <FILES>\n";
		
		if ( $F[6] =~ /\// ) { #pair-end sequence
			$subrun.= "        <FILE checksum=\"$checksum1\" checksum_method=\"$checksum_method\" filename=\"$fastq1\" filetype=\"$filetype\"/>\n";	
			$subrun.= "        <FILE checksum=\"$checksum2\" checksum_method=\"$checksum_method\" filename=\"$fastq2\" filetype=\"$filetype\"/>\n"; #filename --> $F[9]; filetype --> fastq
			}
		else {
			$subrun.= "        <FILE checksum=\"$checksum1\" checksum_method=\"$checksum_method\" filename=\"$fastq1\" filetype=\"$filetype\"/>\n";	
			}
		$subrun .= "      </FILES>\n";
		$subrun .= "    </DATA_BLOCK>\n";
		$subrun .= "    <RUN_ATTRIBUTES>\n";
		$subrun .= "      <RUN_ATTRIBUTE>\n";
		$subrun .= "        <TAG>actual_read_length</TAG>\n"; ##tag --> actual_read_length
		$subrun .= "        <VALUE>$read1_len</VALUE>\n"; ##value --> $F[10]
		$subrun .= "      </RUN_ATTRIBUTE>\n";
		$subrun .= "      <RUN_ATTRIBUTE>\n";
		$subrun .= "        <TAG>run</TAG>\n"; ##tag --> run
		$subrun .= "        <VALUE>$instrument_name</VALUE>\n"; ##value --> $F[3]
		$subrun .= "      <\/RUN_ATTRIBUTE>\n";
		$subrun .= "      <RUN_ATTRIBUTE>\n";	
		$subrun .= "         <TAG>total_bases</TAG>\n"; ##tag --> total_bases
		$subrun .= "         <UNITS>bp</UNITS>\n"; ##units --> bp
		$subrun .= "         <VALUE>$total_bases</VALUE>\n"; ##value --> $F[11]*$F[10]*2;
		$subrun .= "       </RUN_ATTRIBUTE>\n";
		$subrun .= "    </RUN_ATTRIBUTES>\n";
		$subrun .= "  </RUN>\n\n";
		}
	}

my $samplexml=<<SAMPLE;
<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
\n$subsample</SAMPLE_SET>
SAMPLE

my $experimentxml=<<EXPERIMENT;
<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
\n$subexperiment</EXPERIMENT_SET>
EXPERIMENT

my $runxml=<<RUN;
<?xml version="1.0" encoding="UTF-8"?>
<RUN_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
\n$subrun</RUN_SET>
RUN

##output
mkdir "$out_dir/$sra_id.xml";

#open STUDY,">$out_dir/$sra_id.xml/$sra_id.study.xml" or die $!;
#open SAMPLE,">$out_dir/$sra_id.xml/$sra_id.sample.xml" or die $!;
open EXPERIMENT,">$out_dir/$sra_id.xml/$sra_id.experiment.xml" or die $!;
open RUN,">$out_dir/$sra_id.xml/$sra_id.run.xml" or die $!;
open SUBMISSION,">$out_dir/$sra_id.xml/$sra_id.submission.xml" or die $!;
open FASTQ,">$out_dir/$sra_id.fq.txt" or die $!;

#print STUDY $studyxml;
print EXPERIMENT $experimentxml;
print RUN $runxml;
#print SAMPLE $samplexml;
print SUBMISSION $submissionxml;
print FASTQ $fq_list;

close CONFIG;
close FASTQLIST;
#close STUDY;
close EXPERIMENT;
close RUN;
#close SAMPLE;
close SUBMISSION;
close FASTQ;

system( "cd $out_dir && tar -zcf $sra_id.xml.tar.gz $sra_id.xml" ); #pack and compress the xml file

##vlidate xml code
#print STDERR "\nValidating xml code: sometime maybe has connect network error and that's OK!\n1. Validating study.xml code......\n";
#system ( "xmllint --schema $XSD_DIR/SRA.study.xsd $out_dir/$sra_id.xml/$sra_id.study.xml >/dev/null" );

print STDERR "\n1. Validating submission.xml code......\n";
system ( "xmllint --schema $XSD_DIR/SRA.submission.xsd $out_dir/$sra_id.xml/$sra_id.submission.xml >/dev/null" );

#print STDERR "\n3. Validating sample.xml code......\n";
#system ( "xmllint --schema $XSD_DIR/SRA.sample.xsd $out_dir/$sra_id.xml/$sra_id.sample.xml >/dev/null" );

print STDERR "\n2. Validating experiment.xml code......\n";
system ( "xmllint --schema $XSD_DIR/SRA.experiment.xsd $out_dir/$sra_id.xml/$sra_id.experiment.xml >/dev/null" );

print STDERR "\n3. Validating run.xml code......\n";
system ( "xmllint --schema $XSD_DIR/SRA.run.xsd $out_dir/$sra_id.xml/$sra_id.run.xml >/dev/null" );

print STDERR "\nValidating xml code finished! Plaease submit the directory of $sra_id.xml.tar.gz and $sra_id.fq.txt to 3811!!\n";

##sub sample_define
sub sample_define {
my ($sample,$species,$taxon_id) = @_;
my $sample_xml =<<SAMPLE;
  <SAMPLE alias="$sample">
    <SAMPLE_NAME>
      <COMMON_NAME>$species</COMMON_NAME>
      <TAXON_ID>$taxon_id</TAXON_ID>
    </SAMPLE_NAME>
  </SAMPLE>
SAMPLE
return $sample_xml;
}
