#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --account=djones
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --time=2:00:00
#SBATCH --job-name=stringtie_merge

module load system/Miniconda3

source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
conda activate ./condaenv

OUTDIR="output"

srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL \
  gffcompare \
    -T \
    -r "work/genes.gtf" \
    -o "${OUTDIR}/stringtie_merged_curated_gffcompare" \
    "${OUTDIR}/stringtie_merged_curated.gtf"
