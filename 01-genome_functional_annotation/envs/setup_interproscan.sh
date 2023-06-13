#!/usr/bin/env bash

set -euo pipefail

VERSION="5.62-94.0"

mkdir -p work/interproscan
cd work/interproscan
#wget "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/interproscan-${VERSION}-64-bit.tar.gz"
#wget "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/interproscan-${VERSION}-64-bit.tar.gz.md5"

# Recommended checksum to confirm the download was successful:
#md5sum -c interproscan-${VERSION}-64-bit.tar.gz.md5
# Must return *interproscan-5.61-93.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.

#tar -pxvzf "interproscan-${VERSION}-"*-bit.tar.gz

cd "interproscan-${VERSION}"

python3 setup.py -f interproscan.properties
