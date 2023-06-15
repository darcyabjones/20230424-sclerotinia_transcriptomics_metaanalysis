#!/usr/bin/env bash

set -euo pipefail

module load system/Python-3.11.1

python3 -m venv ./gffpal_env

source gffpal_env/bin/activate
python3 -m pip install git+https://github.com/darcyabjones/gffpal.git
