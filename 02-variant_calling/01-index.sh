#!/usr/bin/env bash

set -euo pipefail


mkdir -p work

cp input/Sscl1980-nuclear.fasta work/genome.fasta
bwa index work/genome.fasta
