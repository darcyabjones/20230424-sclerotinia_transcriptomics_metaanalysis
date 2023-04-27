#!/usr/bin/env bash

set -xeuo pipefail

SRA="$1"
REFERENCE="$2"
GTF="$3"
CDS="$4"
HISAT_INDEX="$5"
KALLISTO_INDEX="$6"
OUTDIR="${7:-output}"

mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/"{falco,fastp,strandedness,crams,bigwigs}

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

    rm -rf -- "${TMPPREFIX}/check_strandedness"
    check_strandedness \
        --nreads 500000 \
        --gtf "${GTF}" \
        --transcripts "${CDS}" \
        --kallisto_index "${KALLISTO_INDEX}" \
        --reads_1 "${TMPPREFIX}/${SRA}_trimmed_1.fastq.gz" \
        --reads_2 "${TMPPREFIX}/${SRA}_trimmed_2.fastq.gz" \
        --outdir "${TMPPREFIX}/check_strandedness" \
    > "${TMPPREFIX}/check_strandedness.txt"
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
    rm -rf -- "${TMPPREFIX}/check_strandedness"
    check_strandedness \
        --nreads 500000 \
        --gtf "${GTF}" \
        --transcripts "${CDS}" \
        --kallisto_index "${KALLISTO_INDEX}" \
        --reads_1 "${TMPPREFIX}/${SRA}_trimmed.fastq.gz" \
        --outdir "${TMPPREFIX}/check_strandedness" \
    > "${TMPPREFIX}/check_strandedness.txt" \
    2>&1
fi

# This comes from Rseqc, it'll show up in multiqc
cp "${TMPPREFIX}/check_strandedness/strandedness_check.txt" "${OUTDIR}/strandedness/${SRA}-strandedness_check.txt"

STRANDOUT=$(tail -n 1 "${TMPPREFIX}/check_strandedness.txt")

if [[ "${STRANDOUT}" == "Data is likely FR/fr-secondstrand" ]]
then
    STRAND="FR"
elif [[ "${STRANDOUT}" == 'Data is likely RF/fr-firststrand' ]]
then
    STRAND="RF"
elif [[ "${STRANDOUT}" == 'Data is likely FR/fr-stranded' ]]
then
    STRAND="F"
elif [[ "${STRANDOUT}" == 'Data is likely RF/rf-stranded' ]]
then
    STRAND="R"
elif [[ "${STRANDOUT}" == 'Data is likely unstranded' ]]
then
    STRAND="unstranded ${STRATEGY}"
else
    echo "ERROR: could not determine strandedness" >&2
    cat "${TMPPREFIX}/check_strandedness.txt" >&2
    exit 1
fi

echo "${STRAND}" > "${OUTDIR}/strandedness/${SRA}.txt"

OUTCRAM="${OUTDIR}/crams/${SRA}.cram"
NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"

BNAME="$(basename "${OUTCRAM}")"
mkdir -p "$(dirname "${OUTCRAM}")"

if [ ! -f "${REFERENCE}.fai" ]
then
  samtools faidx "${REFERENCE}"
fi

if [[ "${STRAND}" != "unstranded" ]]
then
    STRAND_PARAM="--rna-strandness '${STRAND}'"
else
    STRAND_PARAM=""
fi

if [ "${STRATEGY}" == "PE" ]
then
    hisat2 \
      -x "${HISAT_INDEX}" \
      -1 "${TMPPREFIX}/${SRA}_trimmed_1.fastq.gz" \
      -2 "${TMPPREFIX}/${SRA}_trimmed_2.fastq.gz" \
      --summary-file "${OUTDIR}/crams/${SRA}-hisat_summary.txt" \
      --downstream-transcriptome-assembly \
      --max-intronlen 20000 \
      --threads "${NCPUS}" \
      --no-unal \
      --rg-id "${SRA}" \
      --rg "SM:${SRA}" \
      ${STRAND_PARAM} \
    2> "${OUTCRAM}.stderr" \
    | samtools view -u --bam --fai-reference "${REFERENCE}.fai" - \
    | samtools fixmate -m - - \
    | samtools sort -@"${NCPUS}" --reference "${REFERENCE}" -O CRAM,embed_ref -o "${OUTCRAM}" -
else
    hisat2 \
      -x "${HISAT_INDEX}" \
      -U "${TMPPREFIX}/${SRA}_trimmed.fastq.gz" \
      --summary-file "${OUTDIR}/crams/${SRA}-hisat_summary.txt" \
      --downstream-transcriptome-assembly \
      --max-intronlen 20000 \
      --threads "${NCPUS}" \
      --no-unal \
      --rg-id "${SRA}" \
      --rg "SM:${SRA}" \
      ${STRAND_PARAM} \
    2> "${OUTCRAM}.stderr" \
    | samtools view -u --bam --fai-reference "${REFERENCE}.fai" - \
    | samtools sort -@"${NCPUS}" --reference "${REFERENCE}" -O CRAM,embed_ref -o "${OUTCRAM}" -
fi

samtools index "${OUTCRAM}"

samtools stats "${OUTCRAM}" --reference "${REFERENCE}" --ref-seq "${REFERENCE}" \
> "${OUTDIR}/crams/${SRA}-samtools_stats.txt"

code/cram_to_bigwig.sh "${REFERENCE}.fai" "${OUTCRAM}" "${OUTDIR}/bigwigs/${SRA}.bw"

if [[ "${STRAND}" != "unstranded" ]]
then
    code/cram_to_bigwig_stranded.sh "${REFERENCE}.fai" "${OUTCRAM}" "${OUTDIR}/bigwigs/${SRA}"
fi
