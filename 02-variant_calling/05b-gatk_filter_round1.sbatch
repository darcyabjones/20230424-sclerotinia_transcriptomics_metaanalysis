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

REFERENCE_GENOME="$(realpath "input/Sscl1980-nuclear.fasta")"
WORKDIR="${PWD}/work/filter_round1"
OUTDIR="${PWD}/output/round1"
mkdir -p "${OUTDIR}" "${WORKDIR}"

bcftools filter \
    --include '(INFO/DP <= 900) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
    -O z \
    -o "${WORKDIR}/SNP-filtered.g.vcf.gz" \
    "${WORKDIR}/SNP.g.vcf.gz"

bcftools index "${WORKDIR}/SNP-filtered.g.vcf.gz"

bcftools filter \
    --include '(INFO/DP <= 900) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
    -O z \
    -o "${WORKDIR}/INDEL-filtered.g.vcf.gz" \
    "${WORKDIR}/INDEL.g.vcf.gz"

bcftools index "${WORKDIR}/INDEL-filtered.g.vcf.gz"

bcftools filter \
    --include '(INFO/DP <= 900) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
    -O z \
    -o "${WORKDIR}/MIXED-filtered.g.vcf.gz" \
    "${WORKDIR}/MIXED.g.vcf.gz"

bcftools index "${WORKDIR}/MIXED-filtered.g.vcf.gz"

#bcftools filter \
#    --include '(INFO/DP <= 900) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
#    -O z \
#    -o "${WORKDIR}/MNP-filtered.g.vcf.gz" \
#    "${WORKDIR}/MNP.g.vcf.gz"
#
#bcftools filter \
#    --include '(INFO/DP <= 700) && (INFO/DP >= 4) && (INFO/MQ >= 40) && (INFO/QD >= 2) && (INFO/FS <= 125) && (INFO/SOR <= 6) && (INFO/MQRankSum >= -12.5) && (INFO/ReadPosRankSum >= -8)' \
#    -O z \
#    -o "${WORKDIR}/SYMBOLIC-filtered.g.vcf.gz" \
#    "${WORKDIR}/SYMBOLIC.g.vcf.gz"


bcftools concat \
    --allow-overlaps \
    -Ou \
    ${WORKDIR}/*-filtered.g.vcf.gz \
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
   -o "${OUTDIR}/filtered.g.vcf.gz" \
   -

bcftools index --tbi "${OUTDIR}/filtered.g.vcf.gz"
