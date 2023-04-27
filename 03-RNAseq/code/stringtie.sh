#!/usr/bin/env bash

set -euo pipefail

SRA="$1"
GTF="$2"
GENOME="$3"
OUTDIR="${4:-output}"

CRAM="${OUTDIR}/crams/${SRA}.cram"
STRAND=$(cat "${OUTDIR}/strandedness/${SRA}.txt")

NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"


if [ "${STRAND}" == "FR" ] || [ "${STRAND}" == "F" ]
then
    STRANDFLAG="--fr"
elif [ "${STRAND}" == "RF" ] || [ "${STRAND}" == "R" ]
then
    STRANDFLAG="--rf"
else
    echo "Not prediting new transcripts as the data are unstranded" >&2
    exit 0
fi

mkdir -p "${OUTDIR}/stringtie_indiv"
stringtie \
    -o "${OUTDIR}/stringtie_indiv/${SRA}.gtf" \
    -p "${NCPUS}" \
    -G "${GTF}" \
    ${STRANDFLAG} \
    --conservative \
    --ref "${GENOME}" \
    "${CRAM}"
