#!/usr/bin/env bash

set -euo pipefail


mkdir -p work

cp input/Sscl1980-nuclear.fasta work/genome.fasta
bwa index work/genome.fasta

awk -F '\t' '$9 ~ /transcript_biotype "mRNA"/ || $3 == "CDS"' input/Sscl1980-mRNA.gtf > work/genes.gtf
