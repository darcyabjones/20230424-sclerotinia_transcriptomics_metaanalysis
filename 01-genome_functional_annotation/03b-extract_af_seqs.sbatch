#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --export=NONE
#SBATCH --time=24:00:00
#SBATCH --job-name=extract

unset PERL5LIB
module load system/Miniconda3
source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
conda activate ./condaenv

set -euo pipefail

find work/alphafold_structures/ -name "*.cif.gz" | code/extract_structure_seq.py --outfile work/all_alphafold_seqs.fasta --kind cif - 
