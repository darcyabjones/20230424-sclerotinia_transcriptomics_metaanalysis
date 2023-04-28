#!/bin/bash --login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --partition=work
#SBATCH --output=%x-%A.stdout
#SBATCH --error=%x-%A.stderr
#SBATCH --export=NONE
#SBATCH --account=y95
#SBATCH --job-name=filter_variants

module load miniconda3/latest
conda activate ./condaenv/

set -euo pipefail

REFERENCE_GENOME="$(realpath "input/ArME14.fasta")"

bcftools filter \
    --include '(INFO/DP <= 900) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
    -O z \
    -o work/gatk_r1_joint_genotypes_SNP-filtered.g.vcf.gz \
    output/gatk_r1_joint_genotypes_SNP.g.vcf.gz

bcftools index work/gatk_r1_joint_genotypes_SNP-filtered.g.vcf.gz

bcftools filter \
    --include '(INFO/DP <= 900) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
    -O z \
    -o work/gatk_r1_joint_genotypes_INDEL-filtered.g.vcf.gz \
    output/gatk_r1_joint_genotypes_INDEL.g.vcf.gz

bcftools index work/gatk_r1_joint_genotypes_INDEL-filtered.g.vcf.gz

bcftools filter \
    --include '(INFO/DP <= 900) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
    -O z \
    -o work/gatk_r1_joint_genotypes_MIXED-filtered.g.vcf.gz \
    output/gatk_r1_joint_genotypes_MIXED.g.vcf.gz

bcftools index work/gatk_r1_joint_genotypes_MIXED-filtered.g.vcf.gz

#bcftools filter \
#    --include '(INFO/DP <= 900) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
#    -O z \
#    -o work/gatk_r1_joint_genotypes_MNP-filtered.g.vcf.gz \
#    output/gatk_r1_joint_genotypes_MNP.g.vcf.gz
#
#bcftools filter \
#    --include '(INFO/DP <= 700) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
#    -O z \
#    -o work/gatk_r1_joint_genotypes_SYMBOLIC-filtered.g.vcf.gz \
#    output/gatk_r1_joint_genotypes_SYMBOLIC.g.vcf.gz


bcftools concat \
    --allow-overlaps \
    -Ou \
    work/gatk_r1_joint_genotypes_*-filtered.g.vcf.gz \
| bcftools norm \
    --fasta-ref "${REFERENCE_GENOME}" \
    --multiallelics +any \
    -Ou \
    - \
| bcftools view \
   --trim-alt-alleles --min-alleles 2 \
   -f .,PASS \
   -Ou \
   - \
| bcftools sort \
   -Oz \
   -o output/gatk_r1_joint_genotypes-filtered.g.vcf.gz \
   -

bcftools index --tbi output/gatk_r1_joint_genotypes-filtered.g.vcf.gz
