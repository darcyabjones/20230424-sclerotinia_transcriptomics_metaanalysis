#!/usr/bin/env bash

set -euo pipefail

SRA="$1"
STRATEGY="$2"
STRAND="$3"
GENOME="$4"
GTF="$5"
OUTDIR="${6:-output}"

CRAM="${OUTDIR}/crams/${SRA}.cram"

NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"


if [ "${STRAND}" == "FR" ] || [ "${STRAND}" == "F" ]
then
    STRANDFLAG="-s 1"
elif [ "${STRAND}" == "RF" ] || [ "${STRAND}" == "R" ]
then
    STRANDFLAG="-s 2"
else
    STRANDFLAG=""
fi

if [ "${STRAND}" == "FR" ] || [ "${STRAND}" == "RF" ] || [ "${STRATEGY}" == "PE" ]
then
    FRAGMENTS="-p --countReadPairs -C"
else
    FRAGMENTS=""
fi

mkdir -p "${OUTDIR}/featurecounts_indiv"

TMPFILE="work/.tmp${SRA}.bam"

trap "rm -f '${TMPFILE}'*" EXIT

MIN_MQ=5

samtools view \
  --reference "${GENOME}" \
  --min-MQ "${MIN_MQ}" \
  -O BAM \
  --write-index \
  -o "${TMPFILE}" \
  "${CRAM}"

featureCounts \
    ${STRANDFLAG} \
    ${FRAGMENTS} \
    -Q "${MIN_MQ}" \
    -t exon \
    -g gene_id \
    -T "${NCPUS}" \
    -a "${GTF}" \
    -o "${OUTDIR}/featurecounts_indiv/${SRA}.txt" \
    "${TMPFILE}"
