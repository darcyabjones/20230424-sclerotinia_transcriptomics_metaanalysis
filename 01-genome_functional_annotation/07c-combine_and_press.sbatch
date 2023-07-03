#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --time=4:00:00
#SBATCH --job-name=index

module load system/Miniconda3

source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
conda activate ./condaenv/

set -euo pipefail

mkdir -p work

cat work/effector_matches_*.stk > work/effector_matches.stk

hmmbuild --amino --hand work/effector_matches.hmm work/effector_matches.stk

sed -i '/^NAME/ s/-i[0-9]$//' work/effector_matches.hmm

hmmpress work/effector_matches.hmm
