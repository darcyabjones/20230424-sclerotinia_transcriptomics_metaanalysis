#!/usr/bin/env bash

set -euo pipefail

SRA="$1"
GENOME="$2"
GTF="$3"
OUTDIR="${4:-output}"

CRAM="${OUTDIR}/crams/${SRA}.cram"
STRAND=$(cat "${OUTDIR}/strandedness/${SRA}.txt")

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

if [ "${STRAND}" == "FR" ] || [ "${STRAND}" == "RF" ] || [ "${STRAND}" == "unstranded PE" ]
then
    FRAGMENTS="-p --countReadPairs -C"
else
    FRAGMENTS=""
fi

mkdir -p "${OUTDIR}/featurecounts_indiv"

TMPDIR="${TMPDIR:-work}/tmp$$"

trap "rm -rf -- '${TMPDIR}'" EXIT

mkdir -p "${TMPDIR}"

samtools view --reference "${GENOME}" -O BAM -o "${TMPDIR}/${SRA}.bam" "${CRAM}"

featureCounts \
    ${STRANDFLAG} \
    ${FRAGMENTS} \
    -t exon \
    -g gene_id \
    -T "${NCPUS}" \
    -a "${GTF}" \
    -o "${OUTDIR}/featurecounts_indiv/${SRA}.txt" \
    "${TMPDIR}/${SRA}.bam"
