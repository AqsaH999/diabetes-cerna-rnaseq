#!/usr/bin/env bash
# 01_download.sh - fetch all raw data + references (run once, in the rnaseq env)
set -euo pipefail
mkdir -p data/raw data/reference data/processed

# 1. Map the 12 chosen GSM samples to SRR run accessions
pysradb gsm-to-srr \
  GSM5009234 GSM5009235 GSM5009236 GSM5009239 GSM5009275 GSM5009280 \
  GSM5009230 GSM5009243 GSM5009250 GSM5009255 GSM5009257 GSM5009267 \
  --saveto data/processed/subset_srr.tsv

# 2. Download FASTQ (single-end SMART-Seq)
cd data/raw
for SRR in $(cut -f2 ../processed/subset_srr.tsv | tail -n +2); do
  echo ">>> $SRR"
  prefetch "$SRR"
  fasterq-dump "$SRR" --split-files --threads 4 --outdir .
  gzip "${SRR}"_*.fastq
  rm -rf "$SRR"
done
cd ../..

# 3. Prebuilt HISAT2 GRCh38 genome index
cd data/reference
curl -L -O https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzf grch38_genome.tar.gz

# 4. Ensembl 110 annotation (matches index chromosome naming)
curl -L -O https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip -f Homo_sapiens.GRCh38.110.gtf.gz

# 5. Gene annotation lookup: gene_id -> gene_name -> biotype
awk -F'\t' '$3=="gene" {
  id=""; name=""; bt="";
  if (match($9, /gene_id "[^"]+"/))      id=substr($9,RSTART+9,RLENGTH-10);
  if (match($9, /gene_name "[^"]+"/))    name=substr($9,RSTART+11,RLENGTH-12);
  if (match($9, /gene_biotype "[^"]+"/)) bt=substr($9,RSTART+14,RLENGTH-15);
  print id"\t"name"\t"bt
}' Homo_sapiens.GRCh38.110.gtf > gene_annotation.tsv
cd ../..

# 6. Published count matrix (Track B)
curl -L -o data/processed/GSE164416_counts.txt.gz \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164416&format=file&file=GSE164416%5FDP%5Fhtseq%5Fcounts%2Etxt%2Egz"
gunzip -f data/processed/GSE164416_counts.txt.gz

# 7. GEO series matrix (ND/T2D metadata)
curl -L -o data/processed/series_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164416/matrix/GSE164416_series_matrix.txt.gz"

echo "Download complete."
