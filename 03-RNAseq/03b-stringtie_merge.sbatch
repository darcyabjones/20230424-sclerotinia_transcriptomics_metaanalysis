#!/bin/bash --login
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --time=2:00:00
#SBATCH --job-name=stringtie_merge

srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL \
  code/stringtie_merge.sh \
    input/Sscl1980-nuclear.fasta \
    input/Sscl1980-mRNA.gtf
