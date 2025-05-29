#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=20000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/Reference/slurm_output/filter-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/Reference/slurm_output/filter-slurm-%j.err

cd /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/Reference/

awk -vRS=">" -vORS="" -vFS="\n" -vOFS="\n" 'NR>1 && $1!~/D02/ {print ">"$0}' /jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_ALL_20170706_transcripts.fasta > IWGSC_v1.1_ALL_20170706_transcripts_noD.fasta

