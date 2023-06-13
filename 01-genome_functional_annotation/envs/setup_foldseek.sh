#!/usr/bin/env bash

set -euo pipefail

# I was getting segfaults with the conda installed foldseek so i switched to a precomiled version.
cd work/

wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz
tar xvzf foldseek-linux-avx2.tar.gz

export PATH=$(pwd)/foldseek/bin/:$PATH
