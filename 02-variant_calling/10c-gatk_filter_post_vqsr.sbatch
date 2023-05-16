#!/bin/bash --login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --partition=workq
#SBATCH --output=%x-%A.stdout
#SBATCH --error=%x-%A.stderr
#SBATCH --export=NONE
#SBATCH --account=djones
#SBATCH --job-name=filter_variants

module load system/Miniconda3
source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"
conda activate ./condaenv/

set -euo pipefail

REFERENCE_GENOME="$(realpath "input/Sscl1980-nuclear.fasta")"
WORKDIR="${PWD}/work/vqsr"
OUTDIR="${PWD}/output/vqsr"
mkdir -p "${OUTDIR}" "${WORKDIR}"

bcftools filter \
    --exclude '(INFO/DP > 4000)' \
    --soft-filter DP_TOO_HIGH \
    -O u \
    "${WORKDIR}/SNP.g.vcf.gz" \
| bcftools filter \
    --exclude '(INFO/DP < 4)' \
    --soft-filter DP_TOO_LOW \
    -O u \
    - \
| bcftools filter \
    --exclude '(INFO/MQ < 50)' \
    --soft-filter MQ_TOO_LOW \
    -O u \
    - \
| bcftools filter \
    --exclude '(INFO/QD < 5)' \
    --soft-filter QD_TOO_LOW \
    -O z \
    -o "${WORKDIR}/SNP-filtered.g.vcf.gz" \
    -

bcftools index "${WORKDIR}/SNP-filtered.g.vcf.gz"

bcftools filter \
    --exclude '(INFO/DP > 3500)' \
    --soft-filter DP_TOO_HIGH \
    -O u \
    "${WORKDIR}/INDEL.g.vcf.gz" \
| bcftools filter \
    --exclude '(INFO/DP < 4)' \
    --soft-filter DP_TOO_LOW \
    -O u \
    - \
| bcftools filter \
    --exclude '(INFO/MQ < 50)' \
    --soft-filter MQ_TOO_LOW \
    -O u \
    - \
| bcftools filter \
    --exclude '(INFO/QD < 5)' \
    --soft-filter QD_TOO_LOW \
    -O z \
    -o "${WORKDIR}/INDEL-filtered.g.vcf.gz" \
    -

bcftools index "${WORKDIR}/INDEL-filtered.g.vcf.gz"

bcftools filter \
    --exclude '(INFO/DP > 3500)' \
    --soft-filter DP_TOO_HIGH \
    -O u \
    "${WORKDIR}/MIXED.g.vcf.gz" \
| bcftools filter \
    --exclude '(INFO/DP < 4)' \
    --soft-filter DP_TOO_LOW \
    -O u \
    - \
| bcftools filter \
    --exclude '(INFO/MQ < 50)' \
    --soft-filter MQ_TOO_LOW \
    -O u \
    - \
| bcftools filter \
    --exclude '(INFO/QD < 5)' \
    --soft-filter QD_TOO_LOW \
    -O z \
    -o "${WORKDIR}/MIXED-filtered.g.vcf.gz" \
    -

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
   -Ou \
   - \
| bcftools sort \
   -Oz \
   -o "${OUTDIR}/final.vcf.gz" \
   -

bcftools index --tbi "${OUTDIR}/final.vcf.gz"
