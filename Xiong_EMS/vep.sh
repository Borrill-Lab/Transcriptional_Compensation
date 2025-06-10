#!/bin/bash
#SBATCH --ntasks 4
#SBATCH -t 00-12:00
#SBATCH -p jic-long
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/Xiong_data/variant_effects/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/Xiong_data/variant_effects/slurm_output/%x.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk

source package /nbi/software/testing/bin/ensembl-vep-91.3
cd /jic/scratch/groups/Philippa-Borrill/Delfi/Xiong_data/

vep --cache --offline --dir_cache "/jic/scratch/groups/Philippa-Borrill/Delfi/.vep" --cache_version 56\
 --species triticum_aestivum -v --canonical --sift b --fork 4 --synonyms "vep_synonyms.txt" --force_overwrite\
 -i "/jic/scratch/groups/Philippa-Borrill/Delfi/Xiong_data/freebayes/merging/all_chromosomes_plus.vcf" -o "variant_effects/variant_effect_output_all_chromosomes.vcf" --vcf