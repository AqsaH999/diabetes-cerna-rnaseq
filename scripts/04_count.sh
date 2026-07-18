#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/tables results/qc/count

# Count reads per gene (single-end -> no -p; -s 0 = unstranded)
featureCounts -T 6 -s 0 \
  -a data/reference/Homo_sapiens.GRCh38.110.gtf \
  -o results/tables/subset_counts.txt \
  data/bam/*.sorted.bam

# summarize assignment stats into a report
multiqc results/tables -o results/qc/count -n multiqc_count
