#!/usr/bin/perl -w

# Script from Philippa Borrill, modified by Marek Glombik and Delfi Dorussen
#
# Aim of script is to run freebayes for multiple samples

#### paths and references:
my $path = '/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1';
my $ref = "$path/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta";
my $input_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/Xiong_data/hisat";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/Xiong_data/freebayes";

### list of chrom to call variants on separately
my $list_of_chrom = "$path/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai";


#############################

#open the list of chromosomes and go through the lines one by one to set off samtools and freebayes
open (INPUT_FILE, "$list_of_chrom") || die "couldn't open the input file $list_of_chrom!";
while (my $line = <INPUT_FILE>) {
chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";
  
#print "\nmy array: @array\n";
print "\narray element 1: @array[0]\n";
  
  
my $chr = $array[0];
  
chdir("$input_dir") or die "couldn't move to specific input directory $input_dir";
  
  
my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel freebayes tasks
#
#SBATCH -p jic-long
#SBATCH -t 12-00:00
#SBATCH -c 4
#SBATCH --mem=100000
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/Xiong_data/freebayes/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/Xiong_data/freebayes/slurm_output/%x.%N.%j.err
SLURM
  
 my $tmp_file = "$output_dir/tmp/freebayes.$chr";
  
  
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $input_dir\n";
  
        print SLURM "set -e\n";
        print SLURM "source package aeee87c4-1923-4732-aca2-f2aff23580cc\n";
        print SLURM "source package 43f9ae9c-a231-41d7-9d07-142ea9637a4d\n";
        print SLURM "samtools merge -R $chr - *.sorted.markdup.bam | freebayes --stdin -f $ref -0 -F 0.87 --min-coverage 10 --use-best-n-alleles 2   > $output_dir/freebayes.$chr.vcf\n";
  
  close SLURM;
  system("sbatch $tmp_file");
 # unlink $tmp_file;
  
}
