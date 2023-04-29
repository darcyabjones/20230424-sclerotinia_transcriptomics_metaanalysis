#!/bin/bash --login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --partition=workq
#SBATCH --output=%x-%A.stdout
#SBATCH --error=%x-%A.stderr
#SBATCH --export=NONE
#SBATCH --account=djones
#SBATCH --job-name=joint_genotyping

module load system/Miniconda3
conda activate ./condaenv/

set -euo pipefail

GVCFS=( work/gatk_haplotype_round1/*.g.vcf.gz )

COMBINED_FILE="/tmp/$$_combined_genotypes.g.vcf.gz"

trap "rm -f ${COMBINED_FILE}" EXIT

srun --export=ALL gatk --java-options "-Xmx7g" CombineGVCFs \
    --reference input/Sscl1980-nuclear.fasta \
    --output "${COMBINED_FILE}" \
    ${GVCFS[@]/#/--variant }


srun --export=ALL gatk --java-options "-Xmx7g" GenotypeGVCFs \
    --reference input/Sscl1980-nuclear.fasta \
    --output work/gatk_combine_round1.g.vcf.gz \
    --call-genotypes true \
    --variant "${COMBINED_FILE}"
