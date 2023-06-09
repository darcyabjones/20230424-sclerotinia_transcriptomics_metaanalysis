set -euo pipefail

if [ $# -eq 0 ] || [ -z "${1:-}" ]
then
  echo "USAGE: $(basename $0) FAI in.cram outfile.bw"
  exit 0
elif [ "${1:-}" == "-" ]
then
  echo "ERROR: Sorry, this tool doesn't support stdin for the FAI" >&2
  exit 1
elif [ "${2:-}" == "-" ]
then
  echo "ERROR: Sorry, this tool doesn't support stdin for the BAM" >&2
  exit 1
elif [ "${3:-}" == "-" ]
then
  echo "ERROR: Sorry, this tool doesn't writing to stdout" >&2
  exit 1
fi

FAI="${1}"
CRAM="${2}"
BW="${3}"

TMPFILE="/tmp/$$-bam_to_bigwig.bedgraph"
trap "rm -f '${TMPFILE}'" EXIT

# Counts the number of mapped primary reads
NREADS=$(samtools view -F 0x104 -c "${CRAM}" | awk '{print $0 / 1000000}')
bedtools genomecov -bga -split -scale "${NREADS}" -ibam "${CRAM}" > "${TMPFILE}"

# on bioconda as ucsc-bedgraphtobigwig
bedGraphToBigWig "${TMPFILE}" "${FAI}" "${BW}"
