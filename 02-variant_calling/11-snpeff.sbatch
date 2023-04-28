#!/bin/bash --login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --partition=work
#SBATCH --output=%x-%A.stdout
#SBATCH --error=%x-%A.stderr
#SBATCH --export=NONE
#SBATCH --account=y95
#SBATCH --job-name=snpeff

module load miniconda3/latest
conda activate ./condaenv/

set -euo pipefail

REFERENCE_GENOME="$(realpath "input/ArME14.fasta")"
REFERENCE_GFF="$(realpath "input/Ascochyta_rabiei_ME14.gff3")"
VCF="output/gatk_r2_joint_genotypes-vqsr_filtered.vcf.gz"

code/snpeff.sh \
    "ArME14" \
    "${REFERENCE_GENOME}" \
    "${REFERENCE_GFF}" \
    "${VCF}"
