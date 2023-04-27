#!/usr/bin/env bash

set -euo pipefail

SRA="$1"
GTF="$2"
GENOME="$3"
OUTDIR="${4:-outdir}"

CRAM="${OUTDIR}/crams/${SRA}.cram"
STRAND=$(cat "${OUTDIR}//strandedness/${SRA}.txt")

NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"


