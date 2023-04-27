#!/usr/bin/env bash

set -euo pipefail

module load system/R-4.1.2_gcc-9.3.0

Rscript code/combine_featurecounts.R output/featurecounts_indiv output/feature_counts.tsv
