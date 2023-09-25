#!/usr/bin/env bash

set -xeuo pipefail

PDBS=$1
CLUSTERS=$2
CLUSTER=$3
OUT=$4

export TMPDIR="work/tmp${CLUSTER}$$"
mkdir -p "${TMPDIR}/pdb"
trap 'rm -rf -- ${TMPDIR}' EXIT

IDS=( $(awk -v COMP="${CLUSTER}" '$3==COMP {print $5}' "${CLUSTERS}" | sort -u) )
if [ "${#IDS[@]}" -lt 2 ]
then
    echo "WARNING: cluster ${CLUSTER} had fewer than 2 members, so i couldn't align it."
    exit 0
fi

for ID in "${IDS[@]}"
do
    cp "${PDBS}/AF-${ID}-"* "${TMPDIR}/pdb"
done

NCPUS=4
export OMP_NUM_THREADS=1
export NUMBA_NUM_THREADS="${NCPUS}"

caretta-cli \
  --threads "${NCPUS}" \
  --output "${TMPDIR}/caretta_results" \
  "${TMPDIR}/pdb"


mkdir -p "${OUT}/${CLUSTER}"
mkdir -p "${OUT}/${CLUSTER}/aligned_pdbs"

sed '/^>/ s/^>AF-\([^-][^-]*\)-.*/>\1/' "${TMPDIR}/caretta_results/result.fasta" > "${OUT}/${CLUSTER}/${CLUSTER}.fasta"

for F in "${TMPDIR}/caretta_results/result_pdb"/*.pdb
do
  cp "${F}" "${OUT}/${CLUSTER}/aligned_pdbs/$(basename ${F/.pdb.pdb/.pdb})" 
done


tar -zcf "${OUT}/${CLUSTER}/aligned_pdbs.tar.gz" "${OUT}/${CLUSTER}/aligned_pdbs"

rm -rf -- "${OUT}/${CLUSTER}/aligned_pdbs"
