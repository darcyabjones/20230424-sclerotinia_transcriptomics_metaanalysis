#!/bin/bash --login

conda activate ./condaenv

set -euo pipefail

if [ $# -lt 4 ]
then
    echo "USAGE: $0 reference outfile variants cram"
    exit 1
fi


REFERENCE="$(realpath "$1")"
shift
OUTFILE="$1"
shift
VARIANTS="$(realpath "$1")"
shift
CRAM="$(realpath "$1")"

mkdir -p "$(dirname "${OUTFILE}")"

gatk BaseRecalibrator \
    --input "${CRAM}" \
    --reference "${REFERENCE}" \
    --known-sites "${VARIANTS}" \
    --output "${OUTFILE}.before.table"


trap "rm -f '${OUTFILE}.bam' '${OUTFILE}.bam.bai'" EXIT

gatk ApplyBQSR \
    --input "${CRAM}" \
    --reference "${REFERENCE}" \
    --bqsr-recal-file "${OUTFILE}.before.table" \
    --output "${OUTFILE}.bam"

gatk BaseRecalibrator \
    --input "${OUTFILE}.bam" \
    --reference "${REFERENCE}" \
    --known-sites "${VARIANTS}" \
    --output "${OUTFILE}.after.table"

gatk AnalyzeCovariates \
    --before-report-file "${OUTFILE}.before.table" \
    --after-report-file "${OUTFILE}.after.table" \
    --plots-report-file "${OUTFILE}.recal.pdf"

samtools view -O CRAM,embed_ref -o "${OUTFILE}.cram" "${OUTFILE}.bam"
samtools index "${OUTFILE}.cram"
