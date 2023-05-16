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
source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"
conda activate ./condaenv/

set -euo pipefail

GVCFS=( work/gatk_haplotype_round1/*.g.vcf.gz )

COMBINED_FILE="${TMPDIR:-work}/$$_combined_genotypes.g.vcf.gz"

trap "rm -f ${COMBINED_FILE}" EXIT

srun --export=ALL gatk --java-options "-Xmx7g" CombineGVCFs \
    --reference work/genome.fasta \
    --output "${COMBINED_FILE}" \
    ${GVCFS[@]/#/--variant }


srun --export=ALL gatk --java-options "-Xmx7g" GenotypeGVCFs \
    --reference work/genome.fasta \
    --output work/gatk_combine_round1.g.vcf.gz \
    --call-genotypes true \
    --variant "${COMBINED_FILE}"
