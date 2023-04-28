#!/usr/bin/env bash

set -euo pipefail

GENOME="$2"
GTF="$3"
OUTDIR="${4:-output}"

NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"

mkdir -p work
ls "${OUTDIR}/stringtie_indiv/"*.gtf > work/stringtie_indiv_files.txt

stringtie \
    --merge \
    -o "${OUTDIR}/stringtie_merged.gtf" \
    -p "${NCPUS}" \
    -G "${GTF}" \
    --ref "${GENOME}" \
    work/stringtie_indiv_files.txt

gffcompare \
    -T \
    -r "${GTF}" \
    -o "${OUTDIR}/stringtie_gffcompare" \
    "${OUTDIR}/stringtie_merged.gtf"
