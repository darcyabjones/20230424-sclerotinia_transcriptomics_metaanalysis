#!/bin/bash --login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --partition=workq
#SBATCH --output=%x-%A.stdout
#SBATCH --error=%x-%A.stderr
#SBATCH --export=NONE
#SBATCH --account=djones
#SBATCH --job-name=snpeff

module load system/Miniconda3
conda activate ./condaenv/

set -euo pipefail

REFERENCE_GENOME="$(realpath "input/Sscl1980-nuclear.fasta")"
REFERENCE_GFF="$(realpath "input/Sscl1980-mRNA.gtf")"
VCF="output/vqsr/filtered.vcf.gz"

code/snpeff.sh \
    "Ss1980" \
    "${REFERENCE_GENOME}" \
    "${REFERENCE_GFF}" \
    "${VCF}"
