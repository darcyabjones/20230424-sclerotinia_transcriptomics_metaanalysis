#!/usr/bin/env bash

set -xeuo pipefail

SRA="$1"
REFERENCE="$2"
OUTDIR="${3:-output}"

mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/"{falco,fastp,crams,bigwigs}

TMPPREFIX="work/.tmp${SRA}"
mkdir -p "${TMPPREFIX}"
trap "rm -rf -- '${TMPPREFIX}'" EXIT

prefetch --progress --resume yes --verify yes --output-file "${TMPPREFIX}/${SRA}.sra" "${SRA}"
fasterq-dump --progress --outdir "${TMPPREFIX}" "${TMPPREFIX}/${SRA}.sra"

if [ -s "${TMPPREFIX}/${SRA}.fastq" ]
then
  STRATEGY=SE
elif [ -s "${TMPPREFIX}/${SRA}_1.fastq" ]
then
  STRATEGY=PE
else
  echo "No fastqs!"
  exit 1
fi

if [ "${STRATEGY}" == "PE" ]
then
    falco --outdir "${OUTDIR}/falco/${SRA}_1" "${TMPPREFIX}/${SRA}_1.fastq"
    falco --outdir "${OUTDIR}/falco/${SRA}_2" "${TMPPREFIX}/${SRA}_2.fastq"
    fastp \
        --in1 "${TMPPREFIX}/${SRA}_1.fastq" \
        --in2 "${TMPPREFIX}/${SRA}_2.fastq" \
        --out1 "${TMPPREFIX}/${SRA}_trimmed_1.fastq.gz" \
        --out2 "${TMPPREFIX}/${SRA}_trimmed_2.fastq.gz" \
        --qualified_quality_phred 5 \
        --length_required 20 \
        --detect_adapter_for_pe \
        --json "${OUTDIR}/fastp/${SRA}-fastp.json" \
        --html "${OUTDIR}/fastp/${SRA}-fastp.html" \
        -R "${SRA}" \
        --thread 3

    rm -f "${TMPPREFIX}/${SRA}_1.fastq" "${TMPPREFIX}/${SRA}_2.fastq"
else
    falco --outdir "${OUTDIR}/falco/${SRA}" "${TMPPREFIX}/${SRA}.fastq"
    fastp \
        --in1 "${TMPPREFIX}/${SRA}.fastq" \
        --out1 "${TMPPREFIX}/${SRA}_trimmed.fastq.gz" \
        --qualified_quality_phred 5 \
        --length_required 20 \
        --json "${OUTDIR}/fastp/${SRA}-fastp.json" \
        --html "${OUTDIR}/fastp/${SRA}-fastp.html" \
        -R "${SRA}" \
        --thread 3

    rm -f "${TMPPREFIX}/${SRA}.fastq"
fi

OUTCRAM="${OUTDIR}/crams/${SRA}.cram"
NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"

BNAME="$(basename "${OUTCRAM}")"
mkdir -p "$(dirname "${OUTCRAM}")"

if [ ! -f "${REFERENCE}.fai" ]
then
  samtools faidx "${REFERENCE}"
fi

RG="@RG\tID:${SRA}\tSM:${SRA}\tPL:ILLUMINA\tPU:${SRA}\tLB:${SRA}"

if [ "${STRATEGY}" == "PE" ]
then
    bwa mem -t ${OMP_NUM_THREADS:-4} -R "${RG}" -M "${REFERENCE}" "${TMPPREFIX}/${SRA}_trimmed_1.fastq.gz" "${TMPPREFIX}/${SRA}_trimmed_2.fastq.gz" \
      | samtools fixmate -u -m - - \
      | samtools sort -u -@2 -T ".samtools_sort$$" - \
      | samtools markdup -@2 --mode s --reference "${REFERENCE}" -O CRAM,embed_ref - "${OUTPREFIX}.cram"

else
    bwa mem -t ${OMP_NUM_THREADS:-4} -R "${RG}" -M "${REFERENCE}" "${TMPPREFIX}/${SRA}_trimmed.fastq.gz" \
      | samtools sort -u -@2 -T ".samtools_sort$$" - \
      | samtools markdup -@2 --mode s --reference "${REFERENCE}" -O CRAM,embed_ref - "${OUTPREFIX}.cram"
fi

samtools index "${OUTCRAM}"

samtools stats "${OUTCRAM}" --reference "${REFERENCE}" --ref-seq "${REFERENCE}" \
> "${OUTDIR}/crams/${SRA}-samtools_stats.txt"

code/cram_to_bigwig.sh "${REFERENCE}.fai" "${OUTCRAM}" "${OUTDIR}/bigwigs/${SRA}.bw"
