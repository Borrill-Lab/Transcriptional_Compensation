#!/bin/bash
#SBATCH -p nbi-long
#SBATCH --ntasks=1
#SBATCH -t 00-12:00
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/variant_effects_filtering/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/variant_effects_filtering/slurm_output/%x.%N.%j.err

cd /jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/variant_effects/

grep synonymous_variant variant_effect_output_all_chromosomes.vcf > /jic/scratch/groups/Philippa-Borrill/Delfi/X201SC22040269-Z01-F001/variant_effects_filtering/variant_effect_output_synonymous.txt
