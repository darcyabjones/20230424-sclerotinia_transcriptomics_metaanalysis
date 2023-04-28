#!/bin/bash

set -euo pipefail

GENOME_NAME=$1
FASTA=$(realpath "$2")
GFF=$(realpath "$3")
VCF=$(realpath "$4")

FASTA_BNAME="$(basename "${FASTA%.*}")"

mkdir -p "work/snpEff/data/${GENOME_NAME}"
echo -e "# ${GENOME_NAME}\n${GENOME_NAME}.genome : ${GENOME_NAME}" > work/snpEff/snpEff.config

cp "${FASTA}" "work/snpEff/data/${GENOME_NAME}/sequences.fa"
awk '!/^#/ && ( ( $3 == "gene" ) || ( $3 == "mRNA" ) || ($3 == "exon" ) || ($3 == "CDS" ) )' "${GFF}" \
| gffread - -T -o "work/snpEff/data/${GENOME_NAME}/genes.gtf" 2> /dev/null

gffread -y "work/snpEff/data/${GENOME_NAME}/protein.fa" -g "work/snpEff/data/${GENOME_NAME}/sequences.fa" --gtf "work/snpEff/data/${GENOME_NAME}/genes.gtf"
gffread -x "work/snpEff/data/${GENOME_NAME}/cds.fa" -g "work/snpEff/data/${GENOME_NAME}/sequences.fa" --gtf "work/snpEff/data/${GENOME_NAME}/genes.gtf"

cd work/snpEff
snpEff build -gtf22 -config ./snpEff.config -dataDir "${PWD}/data" -v "${GENOME_NAME}"  2> /dev/null > /dev/null

bcftools norm \
    --multiallelics - \
    --fasta-ref "${FASTA}" \
    -O v \
    -o in.vcf \
    "${VCF}"

VCF_NOGZ="${VCF%.gz}"
VCF_NOVCF="${VCF_NOGZ%.vcf}"
VCF_NOBCF="${VCF_NOVCF%.bcf}"

snpEff \
  -config snpEff.config \
  -dataDir "data" \
  "${GENOME_NAME}" \
  "in.vcf" \
  > "out.vcf"

bcftools view \
    -O z \
    -o "${VCF_NOBCF}-split-snpeff.vcf.gz" \
    out.vcf

bcftools index --tbi "${VCF_NOBCF}-split-snpeff.vcf.gz"


bcftools view -O v -o in.vcf "${VCF}"

snpEff \
  -config snpEff.config \
  -dataDir "data" \
  "${GENOME_NAME}" \
  "in.vcf" \
  > "out.vcf"

bcftools view \
    -O z \
    -o "${VCF_NOBCF}-snpeff.vcf.gz" \
    out.vcf

bcftools index --tbi "${VCF_NOBCF}-snpeff.vcf.gz"
