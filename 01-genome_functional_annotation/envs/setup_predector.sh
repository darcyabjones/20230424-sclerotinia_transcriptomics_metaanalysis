#!/usr/bin/env bash

set -euo pipefail

#git clone https://github.com/ccdmb/predector envs/predector
#cd envs/predector
#git checkout dev
#git pull --all
#
#cd -

cd envs/sources

ENVIRONMENT="mamba"

unset PERL5LIB

ENVIRONMENT=docker

curl -s "https://raw.githubusercontent.com/ccdmb/predector/1.2.7/install.sh" \
| bash -s "${ENVIRONMENT}" \
    --conda-prefix ../../predector_condaenv \
    -3 signalp-3.0.Linux.tar.Z \
    -4 signalp-4.1g.Linux.tar.gz \
    -5 signalp-5.0b.Linux.tar.gz \
    -6 signalp-6.0g.fast.tar.gz \
    -t targetp-2.0.Linux.tar.gz \
    -d deeploc-1.0.All.tar.gz \
    -m tmhmm-2.0c.Linux.tar.gz \
    -p phobius101_linux.tar.gz

mamba install -c predector -c conda-forge -c bioconda nextflow=22
