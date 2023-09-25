#!/usr/bin/env bash

set -xeuo pipefail

if [ $# -eq 0 ] || [ -z "${1:-}" ]
then
  echo "USAGE: $(basename $0) ALI.fasta SEQS.fasta OUTDIR"
  exit 0
elif [ $# -ne 3 ]
then
  echo "ERROR: need 3 parameters"
  exit 1
fi


ALI=$1
SEQS=$2
OUTDIR=$3

mkdir -p "${OUTDIR}"
unset TMPDIR

TMPDIR="${TMPDIR:-./}tmp$(basename $1)_$(basename $2)"
mkdir -p "${TMPDIR}"

trap "rm -rf -- '${TMPDIR}'" EXIT

NITERS=5

BNAME=$(basename "${ALI%.*}")
ORIG_ALI="${ALI}"

HMM="${TMPDIR}/${BNAME}_orig.hmm"
ORIG_HMM="${HMM}"

hmmbuild --amino --fast --symfrac 0.5 -n "${BNAME}" "${HMM}" "${ALI}"

for i in $(seq 1 "${NITERS}")
do
  STK="${TMPDIR}/out${i}.stk"
  if [ "${i}" -eq "${NITERS}" ]
  then
    OUTPUTS="-A ${STK} --tblout ${OUTDIR}/${BNAME}.tblout --domtblout ${OUTDIR}/${BNAME}.domtblout"
  else
    OUTPUTS="-A ${STK}"
  fi

  hmmsearch \
    -E 0.1 \
    --domE 0.1 \
    --incE 0.00001 \
    --incdomE 0.00001 \
    ${OUTPUTS} \
    --cpu 2 \
    "${HMM}" \
    "${SEQS}" \
  > /dev/null


  (
    cat "${ORIG_ALI}" ;
    esl-reformat fasta "${STK}"
  ) \
  | sed '/^[^>]/s/[-\.[:space:]]//g' \
  > "${STK}.fasta"
  
  ALI="${STK}"
  hmmalign "${ORIG_HMM}" "${STK}.fasta" > "${ALI}"

  if [ "${i}" -ne "${NITERS}" ]
  then
    HMM="${TMPDIR}/out${i}.hmm"
    hmmbuild --amino --hand -n "${BNAME}" "${HMM}" "${ALI}"
  fi
done

cp "${ALI}" "${OUTDIR}/${BNAME}.stk"

