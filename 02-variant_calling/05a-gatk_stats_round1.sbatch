#!/bin/bash --login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=3:00:00
#SBATCH --partition=work
#SBATCH --output=%x-%A.stdout
#SBATCH --error=%x-%A.stderr
#SBATCH --export=NONE
#SBATCH --account=y95
#SBATCH --job-name=variant_stats

module load miniconda3/latest
conda activate ./condaenv/

set -euo pipefail

VARIANT_TYPES=( SNP MNP INDEL SYMBOLIC MIXED )
REFERENCE_GENOME="$(realpath "input/ArME14.fasta")"
REFERENCE_GENES="$(realpath "input/Ascochyta_rabiei_ME14.gff3")"
VARIANTS="$(realpath "work/gatk_r1_joint_genotypes.g.vcf.gz")"

for VARIANT_TYPE in "${VARIANT_TYPES[@]}"
do
    OUTPUT_BNAME="${PWD}/output/$(basename ${VARIANTS%.g.vcf.gz})_${VARIANT_TYPE}"
    mkdir -p "$(basename "${OUTPUT_BNAME}")"
    srun --export=ALL gatk SelectVariants -R "${REFERENCE_GENOME}" --variant "${VARIANTS}" --select-type-to-include "${VARIANT_TYPE}" --output "${OUTPUT_BNAME}.g.vcf.gz"

    NSAMPLES=$(zcat "${OUTPUT_BNAME}.g.vcf.gz" | (grep -v "^#" || :) | wc -l)
    if [ "${NSAMPLES}" -gt 0 ]
    then
        RSCRIPTNAME=/tmp/$$-plot_stats.Rmd
        trap "rm -f ${RSCRIPTNAME}" EXIT
        cp code/plot_stats.Rmd "${RSCRIPTNAME}"
        sed -i "s#INSERT_VARIANT_FILE_HERE_PLEASE#${OUTPUT_BNAME}.g.vcf.gz#" "${RSCRIPTNAME}"
        sed -i "s#INSERT_BIALLELIC_ARG_HERE_PLEASE#FALSE#" "${RSCRIPTNAME}"
        sed -i "s#INSERT_GENOME_FILE_HERE_PLEASE#${REFERENCE_GENOME}#" "${RSCRIPTNAME}"
        sed -i "s#INSERT_GFF_FILE_HERE_PLEASE#${REFERENCE_GENES}#" "${RSCRIPTNAME}"
        sed -i "s#INSERT_GENE_DIST_THRESHOLD_HERE_PLEASE#1000#" "${RSCRIPTNAME}"
        cat "${RSCRIPTNAME}"
        Rscript -e "rmarkdown::render('${RSCRIPTNAME}', output_format='html_document', output_file='${OUTPUT_BNAME}.stats.html', knit_root_dir='${PWD}', clean=TRUE)"
        rm -f "${RSCRIPTNAME}"
    else
        echo "MESSAGE: got no variants for ${VARIANT_TYPE} so skipping stats plotting."
    fi
done
