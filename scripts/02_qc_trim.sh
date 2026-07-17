#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/qc/raw results/qc/trimmed data/trimmed

# 1) QC on raw reads
fastqc data/raw/*_1.fastq.gz -o results/qc/raw -t 4
multiqc results/qc/raw -o results/qc/raw -n multiqc_raw

# 2) Trim each sample with fastp (single-end)
for R1 in data/raw/*_1.fastq.gz; do
  SAMP=$(basename "$R1" _1.fastq.gz)
  echo ">>> trimming $SAMP"
  fastp \
    -i "data/raw/${SAMP}_1.fastq.gz" \
    -o "data/trimmed/${SAMP}_1.trim.fastq.gz" \
    --qualified_quality_phred 20 --length_required 25 \
    --thread 4 \
    --html "results/qc/trimmed/${SAMP}_fastp.html" \
    --json "results/qc/trimmed/${SAMP}_fastp.json"
done

# 3) QC on trimmed reads
fastqc data/trimmed/*_1.trim.fastq.gz -o results/qc/trimmed -t 4
multiqc results/qc/trimmed -o results/qc/trimmed -n multiqc_trimmed
