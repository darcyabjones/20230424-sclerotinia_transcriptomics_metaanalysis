#!/usr/bin/env bash

set -euo pipefail

unset PERL5LIB

mamba create -p ./tmbed_condaenv -c predector -c conda-forge -c bioconda tmbed=1.0.0 
