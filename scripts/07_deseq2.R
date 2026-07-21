# ==========================================================
# 04_deseq2.R  — ND vs T2D differential expression (Track B)
# ==========================================================
library(DESeq2); library(EnhancedVolcano); library(data.table)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# --- 1. Build metadata from the GEO series titles (skip if already made) ---
if (!file.exists("data/processed/metadata.csv")) {
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164416/matrix/GSE164416_series_matrix.txt.gz"
  download.file(url, "data/processed/series_matrix.txt.gz")
  lines  <- readLines(gzfile("data/processed/series_matrix.txt.gz"))
  titles <- gsub('"', '', strsplit(grep("^!Sample_title", lines, value = TRUE), "\t")[[1]][-1])
  meta_all <- data.frame(sample    = sub("Islets_(DP[0-9]+)_.*", "\\1", titles),
                         condition = sub("Islets_DP[0-9]+_(.*)", "\\1", titles))
  write.csv(meta_all, "data/processed/metadata.csv", row.names = FALSE)
}

# --- 2. Load count matrix + metadata ---
cts  <- as.matrix(fread("data/processed/GSE164416_counts.txt"), rownames = 1)
cts  <- cts[grepl("^ENSG", rownames(cts)), ]          # keep real genes only
meta <- read.csv("data/processed/metadata.csv", row.names = 1)

# --- 3. Keep only ND and T2D, match samples to matrix columns ---
meta   <- meta[meta$condition %in% c("ND", "T2D"), , drop = FALSE]
common <- intersect(rownames(meta), colnames(cts))
meta   <- meta[common, , drop = FALSE]
cts    <- cts[, common]
meta$condition <- factor(meta$condition, levels = c("ND", "T2D"))  # ND = reference
mode(cts) <- "integer"
cat("Samples used:\n"); print(table(meta$condition))

# --- 4. Run DESeq2 ---
dds <- DESeqDataSetFromMatrix(cts, meta, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]              # drop near-zero genes
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "T2D", "ND"))
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), "results/tables/deseq2_T2D_vs_ND.csv")

# --- 5. Save a normalized matrix for later (heatmap on Day 7) ---
vsd <- vst(dds, blind = TRUE)
saveRDS(assay(vsd), "results/tables/vst_matrix.rds")

# --- 6. PCA sanity check ---
pdf("results/figures/pca.pdf", width = 6, height = 5)
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

# --- 7. Volcano plot ---
pdf("results/figures/volcano.pdf", width = 7, height = 6)
print(EnhancedVolcano(res, lab = rownames(res),
                      x = "log2FoldChange", y = "padj",
                      pCutoff = 0.05, FCcutoff = 1.0, title = "T2D vs ND"))
dev.off()

# --- 8. Quick summary in the console ---
cat("\nSignificant genes (padj<0.05 & |log2FC|>1):",
    sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm = TRUE), "\n")
summary(res)