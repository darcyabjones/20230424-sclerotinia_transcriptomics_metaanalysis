#!/usr/bin/env bash


set -euo pipefail


curl -OJX GET \
    "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_001857865.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,GENOME_GTF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCA_001857865.1.zip" \
    -H "Accept: application/zip"

unzip -d GCA_001857865.1 GCA_001857865.1.zip

rm -f GCA_001857865.1.zip

rm -rf -- reference_genome
mkdir -p reference_genome

mv GCA_001857865.1/ncbi_dataset/data/GCA_001857865.1/GCA_001857865.1_ASM185786v1_genomic.fna reference_genome/Sscl1980-nuclear.fasta
mv GCA_001857865.1/ncbi_dataset/data/GCA_001857865.1/genomic.gff reference_genome/Sscl1980-mRNA.gff3
mv GCA_001857865.1/ncbi_dataset/data/GCA_001857865.1/genomic.gtf reference_genome/Sscl1980-mRNA.gtf
mv GCA_001857865.1/ncbi_dataset/data/GCA_001857865.1/protein.faa reference_genome/Sscl1980-protein.faa
mv GCA_001857865.1/ncbi_dataset/data/GCA_001857865.1/cds_from_genomic.fna reference_genome/Sscl1980-cds.fna
mv GCA_001857865.1/ncbi_dataset/data/GCA_001857865.1/sequence_report.jsonl reference_genome/
mv GCA_001857865.1/ncbi_dataset/data/assembly_data_report.jsonl reference_genome/

rm -rf -- GCA_001857865.1

awk -v OFS="\t" \
    '{GBK=gensub(/.*genbankAccession":"([^",]+).*/, "\\1", "g", $0); CHR=gensub(/.*chrName":"([^",]+).*/, "\\1", "g", $0); printf "%s\tchr%0.2d\n", GBK, CHR}' \
    reference_genome/sequence_report.jsonl \
> reference_genome/chr_map.tsv

awk -v IF="reference_genome/chr_map.tsv" '
    BEGIN {
        while (getline < IF) {
            split($0,line,"\t")
            IDS[line[1]]=line[2];
        }
    }
    /^>/ {
        GBK=gensub(/>([^[:space:]]+)[[:space:]]Sclerotinia.*/, "\\1", "g", $0);
        CHR=IDS[GBK];
        $0=gensub(/^>(.*)$/, "\\1", "g", $0);
        $0=">"CHR" "$0;
        print
    }
    /^[^>]/ {print}
' reference_genome/Sscl1980-nuclear.fasta \
> reference_genome/Sscl1980-nuclear.fasta.tmp

mv reference_genome/Sscl1980-nuclear.fasta.tmp reference_genome/Sscl1980-nuclear.fasta

bgzip --force reference_genome/Sscl1980-nuclear.fasta
bgzip --force --reindex reference_genome/Sscl1980-nuclear.fasta.gz
samtools faidx reference_genome/Sscl1980-nuclear.fasta.gz

awk -F "\t" -v OFS="\t" -v IF="reference_genome/chr_map.tsv" '
    BEGIN {
        while (getline < IF) {
            split($0,line,"\t")
            IDS[line[1]]=line[2];
        }
    }
    $0 ~ /^##sequence-region/ {
      split($0,line," ");
      sid=IDS[line[2]];
      printf "%s %s %s %s\n", line[1], sid, line[3], line[4];
      next
    }
    $0 ~ /^#/ {print}
    $3 == "region" {next}
    $3 ~ /gene|mRNA|exon|CDS/ {
       $1=IDS[$1];
       $9=gensub(/rna-gnl\|PRJNA348385\|mRNA\./, "", "g", $9);
       print
    }
' reference_genome/Sscl1980-mRNA.gff3 \
> reference_genome/Sscl1980-mRNA.gff3.tmp

(grep '^#' reference_genome/Sscl1980-mRNA.gff3.tmp || :) > reference_genome/Sscl1980-mRNA.gff3.tmp2
(grep -v '^#' reference_genome/Sscl1980-mRNA.gff3.tmp || :) \
 | sort -k1,1 -k4,4n \
 >> reference_genome/Sscl1980-mRNA.gff3.tmp2

mv reference_genome/Sscl1980-mRNA.gff3.tmp2 reference_genome/Sscl1980-mRNA.gff3
rm -f reference_genome/Sscl1980-mRNA.gff3.tmp

bgzip --force reference_genome/Sscl1980-mRNA.gff3
bgzip --force --reindex reference_genome/Sscl1980-mRNA.gff3.gz
tabix -p gff reference_genome/Sscl1980-mRNA.gff3.gz

awk -F "\t" -v OFS="\t" -v IF="reference_genome/chr_map.tsv" '
    BEGIN {
        while (getline < IF) {
            split($0,line,"\t");
            IDS[line[1]]=line[2];
        }
    }
    $0 ~ /^#/ {print}
    $3 == "region" {next}
    $3 ~ /transcript|exon|CDS/ {
       tid=gensub(/^.*orig_transcript_id "gnl\|PRJNA348385\|mRNA\.([^",]+).*$/, "\\1", "g", $9); gid="gene-"tid;
       $9=gensub(/ transcript_id\s+"gnl\|PRJNA348385\|mRNA\.([^"]+)"/, " transcript_id \"\\1\"", "g", $9);
       $9=gensub(/gene_id\s+"[^"]+"/, "gene_id \""gid"\"", "g", $9);
       print $0
    }
' reference_genome/Sscl1980-mRNA.gtf \
> reference_genome/Sscl1980-mRNA.gtf.tmp

(grep '^#' reference_genome/Sscl1980-mRNA.gtf.tmp || :) > reference_genome/Sscl1980-mRNA.gtf.tmp2
(grep -v '^#' reference_genome/Sscl1980-mRNA.gtf.tmp || :) \
 | sort -k1,1 -k4,4n \
 >> reference_genome/Sscl1980-mRNA.gtf.tmp2

mv reference_genome/Sscl1980-mRNA.gtf.tmp2 reference_genome/Sscl1980-mRNA.gtf
rm -f reference_genome/Sscl1980-mRNA.gtf.tmp

bgzip --force reference_genome/Sscl1980-mRNA.gtf
bgzip --force --reindex reference_genome/Sscl1980-mRNA.gtf.gz
tabix -p gff reference_genome/Sscl1980-mRNA.gtf.gz


awk -v IF="reference_genome/chr_map.tsv" '
    BEGIN {
      while (getline < IF) {
        split($0,line,"\t");
        ids[line[1]]=line[2];
      }
    }
    /^>/ {
      $0=gensub(/^[^[:space:]]+\s(.*)$/, "\\1", "g", $0); $0=gensub(/^\[locus_tag=([^\]]+)\]/, ">\\1", "g", $0);
      print
    }
    /^[^>]/ {print}
' reference_genome/Sscl1980-cds.fna \
> reference_genome/Sscl1980-cds.fna.tmp

mv reference_genome/Sscl1980-cds.fna.tmp reference_genome/Sscl1980-cds.fna
bgzip --force reference_genome/Sscl1980-cds.fna
bgzip --force --reindex reference_genome/Sscl1980-cds.fna.gz
samtools faidx reference_genome/Sscl1980-cds.fna.gz

awk '
  /^>/ {
    pid=gensub(/>([^[:space:]]+)[[:space:]].*/, "\\1", "g", $0);
    $0=gensub(/^.*sscle/, ">sscle", "g", $0);
    print $0" [protein_id="pid"]"
  }
  /^[^>]/ {print}
' reference_genome/Sscl1980-protein.faa \
> reference_genome/Sscl1980-protein.faa.tmp

mv reference_genome/Sscl1980-protein.faa.tmp reference_genome/Sscl1980-protein.faa
bgzip --force reference_genome/Sscl1980-protein.faa
bgzip --force --reindex reference_genome/Sscl1980-protein.faa.gz
samtools faidx reference_genome/Sscl1980-protein.faa.gz
