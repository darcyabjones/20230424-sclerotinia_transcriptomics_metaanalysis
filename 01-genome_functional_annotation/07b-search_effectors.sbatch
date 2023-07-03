#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --array=30
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --export=NONE
#SBATCH --time=24:00:00
#SBATCH --job-name=search
#0-29

module load system/Miniconda3

source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
conda activate ./condaenv/

set -euo pipefail

mkdir -p work

INDEX="${SLURM_ARRAY_TASK_ID}"

jackhmmer \
  -o "work/effector_matches_${INDEX}.tsv" \
  -A "work/effector_matches_${INDEX}.stk" \
  --cpu 4 \
  -N 3 \
  "work/effectordb_${INDEX}.fasta" \
  work/uniref50.fasta
