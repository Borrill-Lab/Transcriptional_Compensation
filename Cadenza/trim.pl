#!/usr/bin/perl -w

# Script by Arunkumar Ramesh, modified by Delfi Dorussen
#
# Aim of script is to run trimmomatic for multiple samples

#### paths and references:
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/raw_data/all_reads/";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/trim/";

### lists of samples (text file containing directory/subdirectory with .fastq.gz to trim - text file must be in $input_list_dir):
my $input_list_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/raw_data/all_reads/";
my $input_for_trim = "/jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/raw_data/all_reads/input_for_kallisto.txt";


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$input_list_dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_trim") || die "couldn't open the input file $input_for_trim!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";

#print "\nmy array: @array\n";
#print "\narray element 1: @array[0]\n";

my $sample = $array[0];
my $f1 = $array[1];
my $f2 = $array[2];


chdir("$read_path_triticum") or die "couldn't move to specific read directory $read_path_triticum";


my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel trimmomatic tasks
#
#SBATCH -p nbi-medium
#SBATCH -t 0-04:00
#SBATCH -c 4
#SBATCH --mem=120000
#SBATCH -J trim
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/trim/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/trim/slurm_output/%x.%N.%j.err
SLURM

 my $tmp_file = "$output_dir/tmp/trim.$sample";


  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";

  print SLURM "set -e\n";

  print SLURM "source package 50fcf79b-73a3-4f94-9553-5ed917823423\n";

	print SLURM "trimmomatic PE -phred33 -threads 4 $f1 $f2 $output_dir/$sample"."_1.paired.fq.gz $output_dir/$sample"."_1.unpaired.fq.gz $output_dir/$sample"."_2.paired.fq.gz $output_dir/$sample"."_2.unpaired.fq.gz  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:80\n";

    close SLURM;
  	system("sbatch $tmp_file");
 # unlink $tmp_file;

}

            close(INPUT_FILE);

