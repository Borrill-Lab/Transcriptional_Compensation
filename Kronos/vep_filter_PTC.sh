#!/bin/bash
#SBATCH -p jic-long
#SBATCH --ntasks=1
#SBATCH -t 00-12:00
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/X204SC24093218-Z01-F003/variant_effects_filtering/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/X204SC24093218-Z01-F003/variant_effects_filtering/slurm_output/%x.%N.%j.err

cd /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/X204SC24093218-Z01-F003/variant_effects/

grep stop_gained variant_effect_output_all_chromosomes.vcf > /jic/scratch/groups/Philippa-Borrill/Delfi/Kronos_TA/X204SC24093218-Z01-F003/variant_effects_filtering/variant_effect_output_PTC.txt
