#!/bin/bash --login
#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --export=NONE
#SBATCH --time=4:00:00
#SBATCH --job-name=search

module load system/Miniconda3

source /home/djones/.bashrc && eval "$(conda shell.bash hook)"
conda activate ./condaenv/

set -euo pipefail

mkdir -p work

hmmscan --domtblout work/effector_matches_matches.domtab work/effector_matches.hmm work/Sscl1980.faa

predutils r2js \
  --database-version  2.0 \
  -o output/predector_protein/effectordb_matches.ldjson \
  effectordb work/effector_matches_matches.domtab work/Sscl1980.faa

predutils gff --outfile output/predector_protein/effectordb_matches.gff3 output/predector_protein/effectordb_matches.ldjson 

predutils map_to_genome --split "source" -o "output/predector_genome/effectordb_matches.gff3" work/Sscl1980-mRNA.gff3 output/predector_protein/effectordb_matches.gff3
