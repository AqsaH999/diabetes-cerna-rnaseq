#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/bam results/qc/align
IDX=data/reference/grch38/genome   # HISAT2 index prefix

for R1 in data/trimmed/*_1.trim.fastq.gz; do
  SAMP=$(basename "$R1" _1.trim.fastq.gz)
  echo ">>> aligning $SAMP"
  hisat2 -p 6 -x "$IDX" \
    -U "$R1" \
    --new-summary --summary-file "results/qc/align/${SAMP}.hisat2.log" \
  | samtools sort -@ 4 -o "data/bam/${SAMP}.sorted.bam" -
  samtools index "data/bam/${SAMP}.sorted.bam"
done

# aggregate alignment rates into a report
multiqc results/qc/align -o results/qc/align -n multiqc_align
