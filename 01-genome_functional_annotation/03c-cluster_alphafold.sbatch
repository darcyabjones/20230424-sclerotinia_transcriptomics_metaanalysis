#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --time=8:00:00
#SBATCH --job-name=foldseek

#module load system/Miniconda3

#source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
#conda activate ./condaenv
module load bioinfo/mmseqs2-v14-7e284

set -euo pipefail

export TMPDIR="work/tmp$$"
mkdir -p "${TMPDIR}"
trap 'rm -rf -- ${TMPDIR}' EXIT


export PATH="${PATH:-}:${PWD}/work/foldseek/bin"

"${PWD}"/work/foldseek/bin/foldseek easy-cluster work/alphafold_structures_trimmed work/foldseek_results "${TMPDIR}" \
  -s 9.5 \
  -c 0.7 \
  --cov-mode 0 \
  --cluster-mode 0 \
  -a \
  --threads 16
