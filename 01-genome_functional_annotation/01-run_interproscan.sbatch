#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --export=NONE
#SBATCH --time=10:00:00
#SBATCH --job-name=index

mkdir -p work

module load bioinfo/Java15.0.1

gunzip -c input/Sscl1980-protein.faa.gz > work/Sscl1980.faa

srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL \
  work/interproscan/interproscan-5.62-94.0/interproscan.sh \
    --cpu 4 \
    -f TSV -f XML \
    --output-file-base work/ips \
    --goterms \
    --iprlookup \
    --pathways \
    -i work/Sscl1980.faa
