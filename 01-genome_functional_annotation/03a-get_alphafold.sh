#!/usr/bin/env bash

set -euo pipefail
# Sclerotinia genus
TAXIDS=(
    5179
    5180
    665079
    1095762
    28548
    665077
    1095956
    38451
    665078
    61204
    77105
    1432307
    87234
    90401
    107426
    107579
    262103
    279361
    352851
    504480
    589664
    2036904
    2036905
    2036906
    2075361
    2304652
    2619454
    112198
    126778
    126779
    318626
    318627
    318629
    318870
    364519
    480353
    520545
    660471
    660991
    660992
    728002
    999393
    1137103
    1496367
    1913369
    2831213
)

TAXIDS+=(
    40559
    332648
    999810
    1290391
)


#for TAXID in "${TAXIDS[@]}"
#do
#    echo "Downloading ${TAXID}"
#    gsutil -m cp "gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-${TAXID}-*_v4.tar" . || :
#done


mkdir -p work/alphafold_structures
#mv proteome*.tar work/alphafold_structures

for f in work/alphafold_structures/proteome*.tar
do
    tar --directory=work/alphafold_structures -xf $f
done

find work/alphafold_structures -name "*.tar" -delete
# Have to use this because too many files for rm
find work/alphafold_structures -name "*.json.gz" -delete


