#!/usr/bin/env bash

set -euo pipefail

GENOME="$1"
GTF="$2"
OUTDIR="${3:-output}"

NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"

mkdir -p work
ls "${OUTDIR}/stringtie_indiv/"*.gtf > work/stringtie_indiv_files.txt

# Note that we retain intron retention, as this is much more common in fungi than in other systems.

stringtie \
    --merge \
    -c 5 \
    -F 2 \
    -f 0.05 \
    -i \
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
