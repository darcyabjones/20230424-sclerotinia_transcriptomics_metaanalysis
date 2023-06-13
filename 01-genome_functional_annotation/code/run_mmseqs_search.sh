#!/usr/bin/env bash

set -euo pipefail

QDB=$1
TDB=$2
OUT=$3

module load bioinfo/mmseqs2-v14-7e284

export MDB="work/matches$$/db"
mkdir -p $(dirname "${MDB}")

export TMPDIR="work/tmp$$"
mkdir -p "${TMPDIR}"
trap 'rm -rf -- ${TMPDIR} $(dirname "${MDB}")' EXIT

mkdir -p tmp matches
mmseqs search \
  "${QDB}" \
  "${TDB}" \
  "${MDB}" \
  "tmp" \
  --threads "${NCPUS:-4}" \
  --max-seqs 300 \
  -e 0.01 \
  -s 7 \
  --num-iterations 3 \
  --realign \
  -a

mmseqs convertalis \
  "${QDB}" \
  "${TDB}" \
  "${MDB}" \
  "${OUT}" \
  --threads "${NCPUS:-4}" \
  --format-mode 0 \
  --format-output 'query,target,qstart,qend,qlen,tstart,tend,tlen,evalue,gapopen,pident,alnlen,raw,bits,cigar,mismatch,qcov,tcov'

# predutils r2js \
#   --pipeline-version "${workflow.manifest.version}" \
#   --software-version "${software_version}" \
#   ${db_version_str} \
#   "${database}" search.tsv query/db.fasta \
# > out.ldjson

#rm -rf -- tmp matches search.tsv

