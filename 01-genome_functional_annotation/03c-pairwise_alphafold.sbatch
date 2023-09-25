#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --time=16:00:00
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

"${PWD}"/work/foldseek/bin/foldseek easy-search work/alphafold_structures_trimmed work/alphafold_structures_trimmed work/foldseek_pairwise_results.tsv "${TMPDIR}" \
  -s 9.5 \
  -c 0.8 \
  --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,prob" \
  --cov-mode 2 \
  -a \
  --threads 16
