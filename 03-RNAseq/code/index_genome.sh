#!/usr/bin/env bash

set -euo pipefail

GENOME=$1
CDS=$2
GTF=$3
OUTDIR="${4:-work}"

hisat2_extract_exons.py "${GTF}" > "${OUTDIR}/hisat2_exons.txt"
hisat2_extract_splice_sites.py "${GTF}" > "${OUTDIR}/hisat2_splice_sites.txt"
hisat2-build --exon "${OUTDIR}/hisat2_exons.txt" --ss "${OUTDIR}/hisat2_splice_sites.txt" "${GENOME}" "${OUTDIR}/hisat2_index"
kallisto index -i "${OUTDIR}/kallisto_index" "${CDS}"
