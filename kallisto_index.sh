#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=20000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/Reference/slurm_output/index-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/Reference/slurm_output/index-slurm-%j.err

source package /nbi/software/testing/bin/kallisto-0.46.1

cd /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/Reference/

kallisto index -i IWGSC_v1.1_ALL_20170706_transcripts_noD.fasta_index IWGSC_v1.1_ALL_20170706_transcripts_noD.fasta