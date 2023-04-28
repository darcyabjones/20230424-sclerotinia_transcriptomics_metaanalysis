#!/bin/bash --login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --partition=work
#SBATCH --output=%x-%A.stdout
#SBATCH --error=%x-%A.stderr
#SBATCH --export=NONE
#SBATCH --account=y95
#SBATCH --job-name=joint_genotyping

module load miniconda3/latest
conda activate ./condaenv/

set -euo pipefail

GVCFS=( work/gatk_r2_haplotype_caller/*.g.vcf.gz )

COMBINED_FILE="/tmp/$$_combined_genotypes.g.vcf.gz"

trap "rm -f ${COMBINED_FILE}" EXIT

srun --export=ALL gatk --java-options "-Xmx7g" CombineGVCFs \
    --reference input/ArME14.fasta \
    --output "${COMBINED_FILE}" \
    ${GVCFS[@]/#/--variant }


srun --export=ALL gatk --java-options "-Xmx7g" GenotypeGVCFs \
    --reference input/ArME14.fasta \
    --output work/gatk_r2_joint_genotypes.g.vcf.gz \
    --call-genotypes true \
    --variant "${COMBINED_FILE}"
