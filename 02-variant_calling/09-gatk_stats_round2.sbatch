#!/bin/bash --login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --partition=workq
#SBATCH --output=%x-%A.stdout
#SBATCH --error=%x-%A.stderr
#SBATCH --export=NONE
#SBATCH --account=djones
#SBATCH --job-name=variant_stats

module load system/Miniconda3
source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"
conda activate ./condaenv/

set -euo pipefail

VARIANT_TYPES=( SNP MNP INDEL SYMBOLIC MIXED )
REFERENCE_GENOME="$(realpath "work/genome.fasta")"
REFERENCE_GENES="$(realpath "input/Sscl1980-mRNA.gff3")"
VARIANTS="$(realpath "work/gatk_combine_round2.g.vcf.gz")"

WORKDIR="${PWD}/work/filter_round2"
OUTDIR="${PWD}/output/round2"
mkdir -p "${OUTDIR}" "${WORKDIR}"

for VARIANT_TYPE in "${VARIANT_TYPES[@]}"
do
    OUTPUT_BNAME="${WORKDIR}/${VARIANT_TYPE}"
    mkdir -p "$(basename "${OUTPUT_BNAME}")"
    srun --export=ALL gatk SelectVariants -R "${REFERENCE_GENOME}" --variant "${VARIANTS}" --select-type-to-include "${VARIANT_TYPE}" --output "${OUTPUT_BNAME}.g.vcf.gz"

    NSAMPLES=$(zcat "${OUTPUT_BNAME}.g.vcf.gz" | (grep -v "^#" || :) | wc -l)
    if [ "${NSAMPLES}" -gt 0 ]
    then
        RSCRIPTNAME=/tmp/$$-plot_stats.Rmd
        trap "rm -f ${RSCRIPTNAME}" EXIT
        cp code/plot_gatk_vcf_stats.Rmd "${RSCRIPTNAME}"
        sed -i "s#INSERT_VARIANT_FILE_HERE_PLEASE#${OUTPUT_BNAME}.g.vcf.gz#" "${RSCRIPTNAME}"
        sed -i "s#INSERT_BIALLELIC_ARG_HERE_PLEASE#FALSE#" "${RSCRIPTNAME}"
        sed -i "s#INSERT_GENOME_FILE_HERE_PLEASE#${REFERENCE_GENOME}#" "${RSCRIPTNAME}"
        sed -i "s#INSERT_GFF_FILE_HERE_PLEASE#${REFERENCE_GENES}#" "${RSCRIPTNAME}"
        sed -i "s#INSERT_GENE_DIST_THRESHOLD_HERE_PLEASE#1000#" "${RSCRIPTNAME}"
        #cat "${RSCRIPTNAME}"
        Rscript -e "rmarkdown::render('${RSCRIPTNAME}', output_format='html_document', output_file='${OUTDIR}/${VARIANT_TYPE}-stats.html', knit_root_dir='${PWD}', clean=TRUE)"
        rm -f "${RSCRIPTNAME}"
    else
        echo "MESSAGE: got no variants for ${VARIANT_TYPE} so skipping stats plotting."
    fi
done
