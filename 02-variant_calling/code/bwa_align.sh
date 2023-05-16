#!/bin/bash --login

conda activate ./condaenv

set -euo pipefail

if [ $# -lt 4 ]
then
    echo "USAGE: $0 genome prefix r1 r2"
    exit 1
fi

GENOME="$(realpath "$1")"
SAMPLE="$2"
R1="$(realpath "$3")"
R2="$(realpath "$4")"

mkdir -p "work/${SAMPLE}"
mkdir -p "${PWD}/output/crams"
OUTPREFIX="${PWD}/output/crams/${SAMPLE}"

cd "work/${SAMPLE}"

ID=$(zcat "${R1}" \
  | head -n 1 \
  | sed 's/^.*:\([A-Za-z0-9]*\):\([[:digit:]]*\):[[:digit:]]*:[[:digit:]]*:[[:digit:]]* .*$/\1.\2/') || :
PU="${ID}.${SAMPLE}"
LB="${SAMPLE}"
PL="ILLUMINA"

RG="@RG\tID:${PU}\tSM:${SAMPLE}\tPL:${PL}\tPU:${PU}\tLB:${LB}"

bwa mem -t ${OMP_NUM_THREADS:-4} -R "${RG}" -M genome.fasta "${R1}" "${R2}" \
  | samtools fixmate -u -m - - \
  | samtools sort -u -@2 -T ".samtools_sort$$" - \
  | samtools markdup -@2 --mode s --reference genome.fasta -O CRAM,embed_ref - "${OUTPREFIX}.cram"

samtools index "${OUTPREFIX}.cram"

samtools view -O BAM "${OUTPREFIX}.cram" \
| bedtools genomecov -bga -split -ibam - \
> "genome.bedgraph"

samtools faidx genome.fasta
bedGraphToBigWig genome.bedgraph genome.fasta.fai "${OUTPREFIX}.bw"
