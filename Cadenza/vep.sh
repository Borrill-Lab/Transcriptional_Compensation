#!/bin/bash
#SBATCH --ntasks 4
#SBATCH -t 00-12:00
#SBATCH -p nbi-long
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/variant_effects/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/variant_effects/slurm_output/%x.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk

source package /nbi/software/testing/bin/ensembl-vep-91.3
cd /jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/

vep --cache --offline --dir_cache "/jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/.vep" --cache_version 56\
 --species triticum_aestivum -v --canonical --sift b --fork 4 --synonyms "vep_synonyms.txt" --force_overwrite\
 -i "/jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/freebayes/merging/all_chromosomes_plus.vcf" -o "variant_effects/variant_effect_output_all_chromosomes.vcf" --vcf