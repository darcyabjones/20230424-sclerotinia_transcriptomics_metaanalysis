#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --time=2:00:00
#SBATCH --job-name=multiqc

module load system/Miniconda3

source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
conda activate ./condaenv

mkdir -p output/multiqc

srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL \
  multiqc --outdir output/multiqc output
