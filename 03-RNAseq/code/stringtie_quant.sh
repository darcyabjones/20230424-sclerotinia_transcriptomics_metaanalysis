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
    STRANDFLAG="--fr"
elif [ "${STRAND}" == "RF" ] || [ "${STRAND}" == "R" ]
then
    STRANDFLAG="--rf"
else
    STRANDFLAG=""
fi

mkdir -p "${OUTDIR}/stringtie_quant"
mkdir -p "${OUTDIR}/stringtie_ballgown/${SRA}"
stringtie \
    -e \
    -A "${OUTDIR}/stringtie_quant/${SRA}-gene_abund.tab" \
    -b "${OUTDIR}/stringtie_ballgown/${SRA}" \
    -o "${OUTDIR}/stringtie_quant/${SRA}.gtf" \
    -p "${NCPUS}" \
    -G "${GTF}" \
    ${STRANDFLAG} \
    --ref "${GENOME}" \
    "${CRAM}"
