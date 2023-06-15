#!/usr/bin/env bash

module load system/Python-3.11.1

#cd envs

#tar -zxf sources/deeploc-2.0.All.tar.gz
#wget https://raw.githubusercontent.com/ccdmb/predector/dev/conda/recipes/deeploc2/1.0.0/model.patch
#wget https://raw.githubusercontent.com/ccdmb/predector/dev/conda/recipes/deeploc2/1.0.0/deeploc2.patch
#
#patch deeploc2_package/DeepLoc2/deeploc2.py deeploc2.patch
#patch deeploc2_package/DeepLoc2/model.py model.patch

#cd -

python3 -m venv deeploc2_env
source deeploc2_env/bin/activate

python3 -m pip install wheel
python3 -m pip install ./envs/deeploc2_package
