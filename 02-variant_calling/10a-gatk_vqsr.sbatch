#!/bin/bash --login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --partition=workq
#SBATCH --output=%x-%A.stdout
#SBATCH --error=%x-%A.stderr
#SBATCH --export=NONE
#SBATCH --account=djones
#SBATCH --job-name=filter_variants

module load system/Miniconda3
source /home/djones/.bashrc && unset PERL5LIB && eval "$(conda shell.bash hook)"
conda activate ./condaenv/

set -euo pipefail

REFERENCE_GENOME="$(realpath "work/genome.fasta")"
KNOWN_VCF="output/round1/filtered.g.vcf.gz"
INVCF="work/gatk_combine_round2.g.vcf.gz"


mkdir -p output/vqsr
mkdir -p work/vqsr

gatk MakeSitesOnlyVcf \
    -I "${INVCF}" \
    -O work/vqsr/sites_only.vcf.gz

# The --max-gaussians parameter sets the expected number of clusters in modeling.
# If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
# then the tool will tell you there is insufficient data with a No data found error
# message. In this case, try decrementing the --max-gaussians value.


gatk --java-options "-Xmx7g -Xms7g" VariantRecalibrator \
    -V work/vqsr/sites_only.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -mode INDEL \
    --max-gaussians 4 \
    -resource:bootstap,known=false,training=true,truth=true ${KNOWN_VCF} \
    -O output/vqsr/indels.recal \
    --tranches-file output/vqsr/indels.tranches \
    --rscript-file output/vqsr/indels_plots.R

gatk --java-options "-Xmx5g -Xms5g" \
    ApplyVQSR \
    -V "${INVCF}" \
    --recal-file output/vqsr/indels.recal \
    --tranches-file output/vqsr/indels.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode INDEL \
    -O work/vqsr/indel.recalibrated.vcf.gz


gatk --java-options "-Xmx7g -Xms7g" VariantRecalibrator \
    -V work/vqsr/sites_only.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -mode SNP \
    --max-gaussians 6 \
    -resource:bootstrap,known=false,training=true,truth=true ${KNOWN_VCF} \
    -O output/vqsr/snps.recal \
    --tranches-file output/vqsr/snps.tranches \
    --rscript-file output/vqsr/snps_plots.R

gatk --java-options "-Xmx5g -Xms5g" \
    ApplyVQSR \
    -V work/vqsr/indel.recalibrated.vcf.gz \
    --recal-file output/vqsr/snps.recal \
    --tranches-file output/vqsr/snps.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode SNP \
    -O work/vqsr/snp.recalibrated.vcf.gz


bcftools norm \
    --fasta-ref "${REFERENCE_GENOME}" \
    --multiallelics +any \
    -Ou \
    work/vqsr/snp.recalibrated.vcf.gz \
| bcftools view \
   --trim-alt-alleles --min-alleles 2 \
   -Ou \
   - \
| bcftools sort \
   -Oz \
   -o output/vqsr/filtered.vcf.gz \
   -

bcftools index --tbi output/vqsr/filtered.vcf.gz
