#!/bin/bash --login

conda activate ./condaenv

set -euo pipefail

if [ $# -lt 3 ]
then
    echo "USAGE: $0 reference outfile cram"
    exit 1
fi


REFERENCE="$(realpath "$1")"
shift
OUTFILE="$1"
shift
CRAM="$(realpath "$1")"

mkdir -p "$(dirname "${OUTFILE}")"

gatk HaplotypeCaller \
    --input "${CRAM}" \
    --output "${OUTFILE}" \
    --reference "${REFERENCE}" \
    --sample-ploidy 1 \
    --emit-ref-confidence GVCF

