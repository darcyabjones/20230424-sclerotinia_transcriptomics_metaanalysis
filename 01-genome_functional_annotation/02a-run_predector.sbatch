#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=1-00:00:00
#SBATCH --job-name=predector

mkdir -p work

module load system/Miniconda3

source /home/djones/.bashrc
unset PERL5LIB
eval "$(conda shell.bash hook)"
conda activate ./condaenv

nextflow run \
  -resume \
  -with-conda "${PWD}/predector_condaenv"  \
  -w work/predector_work \
  -ansi-log false \
  envs/predector/main.nf \
    -profile c16,r64 \
    --proteome work/Sscl1980.faa \
    --outdir work/predector_results
  
