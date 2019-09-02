#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd;


main();

sub main{

my ($bam,$threads,$output_dir,$prefix,$path_to_bwa,$samtools_path,$bedtools_path) = process_commandline();
### process input files: fasta and bam format files;
if (@ARGV==0) {
	print "Please provide input files for isomiR seeker(fastq and bam format)! Aborting ...\n";
	exit;
}


my $inputFile=$ARGV[0];

### get current working directory and  scripts directory###
my $current_working_dir = getcwd();
my $scripts_dir = $Bin;

unless ($current_working_dir =~ /\/$/){
      $current_working_dir =~ s/$/\//;
    }
unless ($scripts_dir =~ /\/$/){
      $scripts_dir =~ s/$/\//;
    }

if ($output_dir eq '') {
	
	$output_dir=$current_working_dir;
}

### build isomiR file for output results;

my $isomiR="isomiR";
if(-d "$output_dir$isomiR"){
	deldir("$output_dir$isomiR");
	mkdir "$output_dir$isomiR";
}
else{
	mkdir "$output_dir$isomiR";
}

##################### seeker miR begin #############


chdir ($current_working_dir) or die "Can't move to $current_working_dir: $!";
my $total_miR=0;
if($bam){
#	my $re0=system "rm ${inputFile}.bai";
#	my $re1=system "$samtools_path index $inputFile ${inputFile}.bai";
	my $re=system "$bedtools_path multicov -D -f 1 -r -bed ${scripts_dir}data/hsa_gff3_isomiR.bed -bams $inputFile | cut -f 4,7 | sort -k 1,1 | bedtools groupby -g 1 -c 2 -o sum > $output_dir$isomiR/${prefix}.mature_isomiR.temp.bed"; #>/dev/null 2>&1"; hsa_gff3_isomiR.bed is 0 based; *.bam is 1-based;
	if ($re) {
		exit -1;
	}

my $total_reads=readpipe("samtools view -F 4 -c $inputFile");
### count and RPM ###
chdir ("$output_dir$isomiR") or die "Can't move to $output_dir$isomiR: $!";

open (OUTMIR,'>',"${prefix}.isomiR.count") or die "can not write to the file ${prefix}.isomiR.count: $!";
open (INMIR,"${prefix}.mature_isomiR.temp.bed") or die "can not open the file ${prefix}.mature_isomiR.temp.bed: $!";
open (OUTMIRRPM,'>',"${prefix}.isomiR.rpm") or die "can not write to the file ${prefix}.isomiR.rpm: $!";

while(<INMIR>){
	chomp;
	my @line=split;
#	$total_miR+=$line[1];
	my $rpm=sprintf("%.2f",$line[1]*1000000/$total_reads);
	print OUTMIR "$line[0]\t$prefix\t$line[1]\n";
	print OUTMIRRPM "$line[0]\t$prefix\t$rpm\n";
}
close INMIR;
close OUTMIR;
close OUTMIRRPM;
unlink "${prefix}.mature_isomiR.temp.bed";
}
else{
	die "Fasta format does not support isomiR seeker!";
}

##################### end ##########################


}



sub process_commandline{
	my $help;
	my $version;
	my $bam;
	my $output_dir;
	my $threads;
	my $prefix;
	my $path_to_bwa;
	my $samtools_path;
	my $bedtools_path;

	
	my $command_line = GetOptions('help|man' => \$help,
					'version' => \$version,
					'b|bam' => \$bam,
					'o|output_dir=s' => \$output_dir,
					't|threads=i' => \$threads,
					'p|prefix=s' => \$prefix,
					'path_to_bwa=s' => \$path_to_bwa,
					'samtools_path=s' => \$samtools_path,
					'bedtools_path=s' => \$bedtools_path,	
		

	);

### EXIT ON ERROR if there were errors with any of the supplied options
	 unless ($command_line){
    die "Please respecify command line options\n";
  }

### HELPFILE
	if ($help){
    print_helpfile();
    exit;
  }

	if($version){
	    print "VERSION:0.1.0\n";
		exit;
 }

 ### parse other options
	unless($bam){
	$bam=0;
	}
	else{
	$bam=1;
	}
	
	if($threads){
	if ($threads < 1) {
		die "Number of threads must be a positive integer";
	}
	}
	else{
	$threads=1;
	}
### OUTPUT DIR PATH
  if ($output_dir){
    unless ($output_dir =~ /\/$/){
      $output_dir =~ s/$/\//;
    }
  }
  else{
	 
      $output_dir ='';
    }


### PATH TO bwa
  ### if a special path to bwa was specified we will use that one, otherwise it is assumed that bwa is in the PATH
  if ($path_to_bwa){
	  if($path_to_bwa =~ /bwa$/){
		  if (-e $path_to_bwa){
			# bwa executable found
	}
		  else{
	  die "Could not find an installation of bwa at the location $path_to_bwa. Please respecify\n";
	}

	  }
	  else{
		unless ($path_to_bwa =~ /\/$/){
        $path_to_bwa =~ s/$/\//;
    }
        if (-d $path_to_bwa){
			$path_to_bwa = "${path_to_bwa}bwa";
	
    }
		else{
			die "The path to bwa provided ($path_to_bwa) is invalid (not a directory)!\n";
    }
	  
	  }

  }
  else{
    $path_to_bwa = 'bwa';
    warn "Path to bwa specified as: $path_to_bwa\n"; 
  }

### PATH to samtools

  if (defined $samtools_path){
      # if Samtools was specified as full command
      if ($samtools_path =~ /samtools$/){
		if (-e $samtools_path){
	  # Samtools executable found
	}
		else{
	  die "Could not find an installation of Samtools at the location $samtools_path. Please respecify\n";
	}
      }
      else{
		unless ($samtools_path =~ /\/$/){
			$samtools_path =~ s/$/\//;
	}
		$samtools_path .= 'samtools';
   		if (-e $samtools_path){
			# Samtools executable found
	}
		else{
	  die "Could not find an installation of Samtools at the location $samtools_path. Please respecify\n";
	}
      }

      warn "Samtools path provided as: '$samtools_path'\n";
    }
    # Check whether Samtools is in the PATH if no path was supplied by the user
  else{
      if (!system "which samtools >/dev/null 2>&1"){ # STDOUT is binned, STDERR is redirected to STDOUT. Returns 0 if samtools is in the PATH
		$samtools_path = `which samtools`;
		chomp $samtools_path;
		warn "Samtools found here: '$samtools_path'\n";
      }
    }

  unless (defined $samtools_path){
      warn "Did not find Samtools on the system.\n";
	  exit;
    }
    
### PATH to bedtools
  if (defined $bedtools_path){
      # if bedtools was specified as full command
      if ($bedtools_path =~ /bedtools$/){
		if (-e $bedtools_path){
	  # Samtools executable found
	}
		else{
	  die "Could not find an installation of bedtools at the location $bedtools_path. Please respecify\n";
	}
      }
      else{
		unless ($bedtools_path =~ /\/$/){
			$bedtools_path =~ s/$/\//;
	}
		$bedtools_path .= 'bedtools';
   		if (-e $bedtools_path){
			# bedtools executable found
	}
		else{
	  die "Could not find an installation of bedtools at the location $bedtools_path. Please respecify\n";
	}
      }

      warn "bedtools path provided as: '$bedtools_path'\n";
    }
    # Check whether bedtools is in the PATH if no path was supplied by the user
  else{
      if (!system "which bedtools >/dev/null 2>&1"){ # STDOUT is binned, STDERR is redirected to STDOUT. Returns 0 if bedtools is in the PATH
		$bedtools_path = `which bedtools`;
		chomp $bedtools_path;
		warn "bedtools found here: '$bedtools_path'\n";
      }
    }

  unless (defined $bedtools_path){
      warn "Did not find bedtools on the system.\n";
	  exit;
    }


## PREFIX FOR OUTPUT FILES
  if ($prefix){
    # removing trailing dots

    $prefix =~ s/\.+$//;

    warn "Using the following prefix for output files: $prefix\n\n";
    sleep(1);
  }
  else{
  $prefix="isomiRseeker";
  }


return($bam,$threads,$output_dir,$prefix,$path_to_bwa,$samtools_path,$bedtools_path);
}


sub deldir{ 

	my ($del_dir) = $_[0]; 
	my (@direct); 
	my (@files); 
	opendir (DIR2,"$del_dir"); 
	my (@allfile)=readdir(DIR2); 
	close (DIR2); 
	foreach (@allfile){ 
		if (-d "$del_dir/$_"){ 
			push(@direct,"$_"); 
} 
		else { 
			push(@files,"$_"); 
} 
} 

my	$files=@files; 
my	$direct=@direct; 
	if ($files ne "0"){ 
		foreach (@files){ 
		unlink("$del_dir/$_"); 
} 
} 
	if ($direct ne "0"){ 
		foreach (@direct){ 
		&deldir("$del_dir/$_") if($_ ne "." && $_ ne ".."); 
} 
} 

	rmdir ("$del_dir"); 
}

