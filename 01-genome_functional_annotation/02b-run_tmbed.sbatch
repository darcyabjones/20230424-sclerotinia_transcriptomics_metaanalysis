#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --time=24:00:00
#SBATCH --job-name=tmbed

module load system/Miniconda3

source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
conda activate ./tmbed_condaenv/

set -euo pipefail

export TMPDIR="work/tmp$$"
mkdir -p "${TMPDIR}"
trap 'rm -rf -- ${TMPDIR}' EXIT

#tmbed embed -f work/Sscl1980.faa -e output/t5_embeddings.h5
tmbed predict -f work/Sscl1980.faa -e output/t5_embeddings.h5 -p output/tmbed.3line --out-format 0
