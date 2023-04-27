#!/bin/bash --login
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --time=2:00:00
#SBATCH --job-name=index

srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL \
  code/index_genome.sh \
    input/Sscl1980-nuclear.fasta \
    input/Sscl1980-cds.fna \
    input/Sscl1980-mRNA.gtf