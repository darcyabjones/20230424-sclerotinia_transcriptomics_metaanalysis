#!/usr/bin/env bash

set -euo pipefail

SRA="$1"
shift
STRAND="$1"
shift
GENOME="$1"
shift
GTF="$1"
shift
CRAMS=( "$@" )
OUTDIR="output"


NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"


if [ "${STRAND}" == "FR" ] || [ "${STRAND}" == "F" ]
then
    STRANDFLAG="--fr"
elif [ "${STRAND}" == "RF" ] || [ "${STRAND}" == "R" ]
then
    STRANDFLAG="--rf"
else
    echo "ERROR: got unknown strand $STRAND" >&2
    exit 1
fi


CRAM="work/.tmp${SRA}.cram"
trap "rm -rf -- '${CRAM}'*" EXIT

samtools merge \
  -f \
  -u \
  --threads "${NCPUS}" \
  --reference "${GENOME}" \
  -o - \
  "${CRAMS[@]}" \
| samtools view \
  --reference "${GENOME}" \
  --min-MQ 5 \
  -O CRAM,embed_ref \
  --write-index -o "${CRAM}" \
  -

mkdir -p "${OUTDIR}/stringtie_indiv"
stringtie \
    -o "${OUTDIR}/stringtie_indiv/${SRA}.gtf" \
    -p "${NCPUS}" \
    -G "${GTF}" \
    ${STRANDFLAG} \
    -g 0 \
    -f 0.1 \
    -j 6 \
    -c 4 \
    -s 4 \
    --ref "${GENOME}" \
    "${CRAM}"
