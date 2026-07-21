#!/usr/bin/env bash
# run_all.sh - full pipeline, in order.
# Shell/Python steps: run with the 'rnaseq' conda env active.
# R steps: require R + Bioconductor packages (see docs/methods.md). Rscript must
#          point to that R install; alternatively run each R script in RStudio
#          after opening diabetes-cerna-rnaseq.Rproj.
set -euo pipefail
cd "$(dirname "$0")"

echo "== Track A: raw data -> validated counts =="
bash   scripts/01_download.sh
bash   scripts/02_qc_trim.sh
bash   scripts/03_align.sh
bash   scripts/04_count.sh
python scripts/05_validate.py
python scripts/06_scatter.py

echo "== Track B: differential expression -> ceRNA network =="
Rscript scripts/07_deseq2.R
Rscript scripts/08_figures.R
Rscript scripts/09_biotype_split.R
Rscript scripts/10_enrichment.R
Rscript scripts/11_cerna_predict.R
Rscript scripts/12_network.R
Rscript scripts/13_network_viz.R
Rscript scripts/14_cerna_correlation.R
Rscript scripts/15_network_supported.R

echo "Pipeline complete. See results/ for tables and figures."
